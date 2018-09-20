
function MID_SSVEP (SubjCode, runNumber)

% Inputs: 'SubjCode', anonymised unique subject code as a string, eg. 'RM'
%         runNumber, scalar giving the current run number.
% Displays 4 interleaved conditions at random intervals (SOAs):
% 1. CD towards
% 2. CD away
% 3. IOVD towards
% 4. IOVD away
% SOAs = stimulus onset asynchronies, the time between the ONSETS of the probes.
% Subjects respond 'up' or '8' for away and 'down' or '2' for towards.
% Runs currently set to 5 minutes (RunDurationSec=300).

% towards: negative sawtooth wave
% away: positive sawtooth wave.


%%
% Set default values if not parsed.
if nargin < 2
    runNumber = 1;
end

if nargin < 1
    SubjCode = 'test';
end

%%
%%%%-------------------------------%%%%
%           Define parameters:
%%%%-------------------------------%%%%

% If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = false;

% Define some of the display parameters:
PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)

% Switch for whether you want the annulus superimposed over the dots:
DrawAnnulus = true;

% Draw the fixation lock/rings
DrawRings = true;

% Spatiotemporal parameters for IOVD and CD:
% Note that these must be precisely the same as those in the pre-generated stimulus files.
% This is more to define the appropriate file names for loading, since the frequency itself will be stored in the parameters struct.
parameters.IOVD_freq = 4; %in Hz
parameters.CD_freq = 4;

% Define the max. change in disparity (horiz shift of dots) for the 2 cues:
% the numerator is in arcmin (but converted to pixels) akin to the amplitude of the sine wave
parameters.IOVD_disparityPixels = 128/60 * PPD;
parameters.CD_disparityPixels = 16/60 * PPD;
% Store the values in ArcMin.
parameters.IOVD_disparityArcMin = 128;
parameters.CD_disparityArcMin = 16;
%%
% work out output filenames:
dateString = datestr(now);
dateString(dateString == ' ') =  '_';
dateString(dateString == '-') =  '_';
dateString(dateString == ':') =  '_';
dateString = dateString(1:end-9); % get rid of the timestamp off the end, don't really need that.
%Set aside information
subjectData.experimentDescriptor = 'MID_SSVEP';
subjectData.subjectCode = SubjCode;
subjectData.runNumber = runNumber;

% File name for the data to be saved as:
subjectData.filename = fullfile('data', [SubjCode, '_', ...  %maybe put the /data prefix in later.
    subjectData.experimentDescriptor, ...
    ['_' , dateString, '_'] ...
    num2str( runNumber ), ...
    '.mat']);

% Just abort if the file already exists:
if exist(subjectData.filename,'file')
    userResponse = input('WARNING: Run file exists. Overwrite? Enter y or n: ','s');
    if ~strcmp( userResponse, 'y' )
        subjectData = [];
        error('Aborting function! File exists!');
    end
end

if ~ismac
    jheapcl; % clear the java heap space.
end

% Choose the screen: it is usually the max screen no. available.
% Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
% So if we're on that machine, need to -1 from the screen number:
[~, CompName] = system('hostname');                                 % find out the computer name
if strncmpi(CompName, 'pspclabrm01', length('pspclabrm01')) ... % normal strcmp not working here, can't figure out why...
        && (length(Screen('Screens'))>1)                            % and there is more than 1 display connected...
    WhichScreen = max( Screen( 'Screens' ) )-1;
else
    WhichScreen = max( Screen( 'Screens' ) );                       % should be right for any other machine!
end

screenRect = Screen('Rect',WhichScreen);                            % get the screen resolution.
centreX = screenRect(3)/2;
centreY = screenRect(4)/2;
RefreshRate = Screen('NominalFrameRate', WhichScreen);

% the effective frame rate PER EYE: since each eye is stimulated on successive video frames
% Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
PerEyeFR = RefreshRate/2; %Per eye f.r. is always half the absolute f.r.

% Save the information about the PC hardware:
subjectData.HostPC = CompName;

%%%%-------------------------------%%%%
%       define response keys
%%%%-------------------------------%%%%

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

% Subjects can press '8' or 'up' for 'away' motion on the numeric keypad.
AwayMotion = [KbName('UpArrow'), KbName('8')];
% Or press '2' or 'down' for 'towards' motion.
TowardsMotion = [KbName('DownArrow'), KbName('2')];
% Finally, press 'q' at any time to quit the program.
RespQuit = KbName('q'); % q to quit.
% NOTE: if you want to change these to the keyboard arrows, the codes are:
% KbName('UpArrow') and KbName('DownArrow').
% Incidentally, these are the same codes as 8 and 2 on the keypad if Num Lock is off.

% Set up a 'shortcut' list of values to be scanned by KbCheck.
% These are the only keyboard responses we are interested in.
% This should save a little bit of time with each call to 'KbCheck', because it only needs to scan these values, and ignore others
scanListVals = [TowardsMotion, AwayMotion, RespQuit]; % 8 for away motion, 2 for towards, q to quit.
scanList = zeros(1,256); scanList(scanListVals) = 1;

%%
%%%%----------------------------------------%%%%
%           Define timing parameters:
%%%%----------------------------------------%%%%

%%%% --------------------------------------------------------%%%%
% *** A note on defining the timing and rounding errors ***
%%%% --------------------------------------------------------%%%%
% We get rounding errors later on when matching the SOAs with the single (per-eye) noise frames (to add the disparity)!
% This is because the duration of 1 per-eye frame (or single-eye frame) has a repeating decimal expansion eg. 0.01666667....
% We could maybe avoid this by defining everything in FRAMES here from the outset instead of SEC
% BUT it is more intuitive for later (analysis, etc) if things are first defined in seconds.
% Instead, we will first round to the nearest second, and then add some jitter to that within a pre-defined range, given below.
% Later, when the disparity is added, we convert the time in sec to integer frames.
% These rounding errors were very strange and I couldn't figure them out, on occasions matlab would say, eg (0.33333==0.33333) = 0.
% Maybe had something to do with 1 value being off somewhere in the 64-bit double-precision arrays... (RM: 17/1/17).

% The possible ONSET times of per-eye frame doublets from 0:1 second.
%PerEyeFrameOnsets1Sec = 0:(1/PerEyeFR):1-(1/PerEyeFR);

% Determine the stimulus ONSETS of each of the probes.
ProbeLengthSec = 0.25; % Duration of the MID probes in sec.
RunDurationSec = 300; % five minute run seems about right for EEG.

% Work out the SOAs of MID target events, drawn at random from a uniform distribution
LB = 1;  % Lower bound on distribution (0.25 s stimu duration + minimum desired gap)
UB = 3;    % Upper bound on distribution

% Establish the stimulus onset asynchronies (SOAs) for each of the MID probes.
cSOAs = []; % the cumulative SOAs
SOAcounter = 1; % Start at 1 because first SOA is determined outside the loop below
% Generate the very first target onset time:
X = LB + (UB - LB) .* rand; %This gives us the actual time in seconds.

cSOAs(1) = X; % cSOAs contains the cumulative SOAs across the whole run for the probes
while cSOAs(SOAcounter) < RunDurationSec-ProbeLengthSec
    SOAcounter = SOAcounter+1;
    % Simulate uniform draws from the interval LB to UB
    X = LB + (UB - LB) .* rand;
    cSOAs(SOAcounter) = cSOAs(SOAcounter-1) + X;
end

% Get rid of any values that might overlap with the end of the run:
cSOAs = cSOAs(cSOAs <= RunDurationSec-ProbeLengthSec);

% Let's assign each of the SOAs to one of the 4 conditions now. We want each to occur with equal probability, so
% let's cull cSOAs back until it's an exact multiple of 4:
while mod(length(cSOAs), 4) > 0
    cSOAs = cSOAs(1:end-1);
end
%figure; bar(diff(cSOAs)); % plot distr of values before culling

% Now we can set up a vector that assigns each one of the cSOAs to one of the 4 conditions,
% in random order with equal probability
% The conditions are coded thus:
% 1 = CD towards
% 2 = CD away
% 3 IOVD towards
% 4 IOVD away
AssignCond = repmat(1:4, 1, length(cSOAs)/4);
AssignCond = AssignCond(randperm(length(AssignCond)));
% Also insert 5 for the blank periods:
AssignCond = [AssignCond; ones(1,length(AssignCond))*5];
AssignCond = reshape(AssignCond, 1, length(AssignCond)*2);
AssignCond = [5 AssignCond]; % and insert a 5 for the start blank

% Now we need to insert the end times of each stimulus (in other words, the ON times of each blank period).
cSOAs(2,:) = cSOAs+ProbeLengthSec;
% Reshape the matrix into a single vector of times:
cSOAs = reshape(cSOAs, 1, length(cSOAs)*2);
% Insert start and end points of the run:
cSOAs = [0 cSOAs RunDurationSec];

% Now, generate a vector of the DURATIONS of each probe (stimulus or blank)
% This can be used to control the length of each blank or probe period.
stimDurations = diff(cSOAs);
% Note that sum(stimDurations) = RunDurationSec

% Set up matrix to store button press responses, and indicate their accuracy.
% We will only check for responses after each pair of PER EYE frames, so the
% length of this matrix will be the equivalent of the total per eye frames.
% Col 1 will store the type of stimulus, col 2 will code a button press as correct or incorrect.
% The type of stimulus will be recorded with every probe, even if no button press is made.
Responses = nan(RunDurationSec*PerEyeFR, 2);

% set aside for saving:
parameters.RunDurationSec = RunDurationSec;
parameters.ProbeLengthSec = ProbeLengthSec;
parameters.ConditionOrder = AssignCond;
parameters.PPD = PPD;
parameters.cSOAs = cSOAs; %including blanks now
parameters.PerEyeFR = PerEyeFR;
% This value determines the timing of the entire experiment: should be 18000 per-eye frames for a 5-min run
parameters.TotalPerEyeFrames = RunDurationSec*PerEyeFR;
parameters.stimDurations = stimDurations; % All stimulus durations, including the blanks.

%%
%%%%-------------------------------%%%%
%           Define stimuli:
%%%%-------------------------------%%%%

% Note that we are now no longer using temporal cosine ramps for the stimuli,
% so we are no longer defining the peak dot contrast (this will be 1 by default).

% Define the dot texture, a square-shaped sheet of dots.
% Make the texture the same size as the height of the screen
% (well really we are using many little dot textures, but 'imsize' defines the region they inhabit)
% Note, this should really be the same for both CD and IOVD...
imsize = screenRect(4);

% specify dot size:
dot_sigma_in_degrees = 0.05;            % size of SD (sigma) of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; % sigma in pixels
dotsize = round(dot_sigma * 10);        % make the dots some multiple of sigma
% NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
% It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.

% make the Gaussian dot profile:
x = (1-dotsize)/2:(dotsize-1)/2;
%x = x/PPD; %rescale x into vis angle
[x,y] = meshgrid(x,x);
y = -y;
[~,r] = cart2pol(x,y);

% This gives us a white-peaked dot (+ve contrast polarity)
env = exp(-r.^2/(2*dot_sigma.^2)); % gaussian window
env = env./max(max(abs(env)));     % normalize peak to +/- 1
env2 = -env;                       % make a negative version to give us a black dot (-ve contrast polarity)

% set up the raised cosine annular window. We need this to be at least as big as the screen,
% so that the texture doesn't poke out from behind it.
% specify parameters for the annulus:
inrad = PPD * 1;     % inner radius of annulus (result in pixels), for fixation spot
outrad = PPD * 6;     % outer radius of annulus (result in pixels)
% define extent of spatial raised cosine at edge of aperture (in pixels)
cos_smooth = dotsize;   % make it one dot size wide
% This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2;     % double the texture size
x0 = (imsize2+1)/2;
y0 = (imsize2+1)/2;
J = ones(imsize2);
for (ii=1:imsize2)
    for (jj=1:imsize2)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > outrad^2)
            J(ii,jj) = 0;
        elseif (r2 < inrad^2)
            J(ii,jj) = 0;
        elseif (r2 > (outrad - cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
        elseif (r2 < (inrad + cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
        end
    end
end

%Set aside some more of the stimulus parameters:
parameters.DotSigmaInDeg = dot_sigma_in_degrees;
parameters.imsize = imsize;
parameters.PeakDotContrast = 1; %(by default)
parameters.DotSigmaInDeg = dot_sigma_in_degrees;

%%
%%%%-------------------------------%%%%
%       Set up the fixation cross:
%%%%-------------------------------%%%%

% Set up the fixation cross or spot:
% This is drawn directly to the screen using Screen('FillRect')
% if you're using a cross instead:
crossWidth = 2;
crossHeight = 10;
fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

% Make the fixation lock rings:
% We have an inner one around fixation and an outer one right on the edge of screen.
% These could probably be defined as a single texture (rather than 2) but I think that will complicate matters with the alpha-blending settings.
% (they are complicated enough already)
ringRadiusInner = PPD*0.4;                  % ring surrounding fixation
ringRadiusOuter = screenRect(4)/2;          % outer edge (radius) of the ring: the edge of the screen
ringWidthInner = ringRadiusInner - PPD/5;   % 1/5 of a degree thick
ringWidthOuter = ringRadiusOuter - PPD/3;   % 1/3 of a degree thick

% Make the rings. Both are in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
% Define the ring:
xx = (1-imsize)/2:(imsize-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
% make the alpha mask for the rings, inner and outer.
ring_alphaInner = ((r>ringWidthInner+1) & (r<ringRadiusInner-1)); % Make the alpha mask a tiny bit thinner than the ring itself.
ring_alphaOuter = ((r>ringWidthOuter+1) & (r<ringRadiusOuter-1));

%%
%%%%-------------------------------%%%%
%           Set up the screen:
%%%%-------------------------------%%%%
ForcedQuit = false; % this is a flag for the exit function to indicate whether the program was aborted
try % Start a try/catch statement, in case something goes awry with the PTB functions
    
    % initialization of the display
    AssertOpenGL;
    % Open PTB onscreen window: We request a 32 bit per colour component
    % floating point framebuffer if it supports alpha-blending. Otherwise
    % the system shall fall back to a 16 bit per colour component framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % required for gamma correction through the PsychImaging pipeline:
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
    % Set the color range to be normalised between 0 and 1 (rather than 0-255):
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
    
    % Open an on screen (grey) window and configure the imaging pipeline
    % Info about the 'blueline' mechanism for synching to the 3D glasses:
    % There seems to be a blueline generation bug on some OpenGL systems.
    % SetStereoBlueLineSyncParameters(windowPtr, windowRect(4)) corrects the
    % bug on some systems, but breaks on other systems.
    % We'll just disable automatic blueline, and manually draw our own bluelines!
    if useHardwareStereo
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5, [], [], [], 1); %flag of 1 engages stereomode
        SetStereoBlueLineSyncParameters(win, windowRect(4)+10);
    else
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);
    end
    
    % Initialise the Vpixx device:
    if UsingVP                      % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        % The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        Datapixx('EnableVideoStereoBlueline');
        Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses
        
        if Datapixx('IsViewpixx3D')     % If it's the Viewpixx3D
            
            % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
            % SHOULD NOT use Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
            % The PROpixx device should have a linear lUT built in, but we will add this here for completeness.
            % The values are the K channel (all RGB) gamma function fits as measured through the goggles, 2/10/15.
            % The spectral measurements were made using the Jaz spectrometer on the ShuttlePC
            K_gamma_left = 2.9712;
            K_gamma_right = 3.4363;
            % We'll use the average of the right and left gamma values
            PsychColorCorrection('SetEncodingGamma', win, 1/mean([K_gamma_left, K_gamma_right]));
            
            %Datapixx('EnableVideoLcd3D60Hz');
            Datapixx('DisableVideoLcd3D60Hz'); %=> weirdly, this seems to cause less crosstalk (see DBaker)
            parameters.DisplayType = 'Viewpixx3D'; % set aside the device type for reference
            Datapixx('RegWr');
            
        elseif Datapixx('IsPropixx')    % if it's the Propixx DLP projector
            
            parameters.DisplayType = 'PROpixx';         % set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram',0); % set to normal RGB video processing for driving the LEDs & DLP MMDs
            % Datapixx('RegWr');
            Datapixx('RegWrRd'); % seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        end
    end
    % No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
    %Define the 'blue line' parameters
    blueRectLeftOn   = [0,                 windowRect(4)-1, windowRect(3)/4,   windowRect(4)];
    blueRectLeftOff  = [windowRect(3)/4,   windowRect(4)-1, windowRect(3),     windowRect(4)];
    blueRectRightOn  = [0,                 windowRect(4)-1, windowRect(3)*3/4, windowRect(4)];
    blueRectRightOff = [windowRect(3)*3/4, windowRect(4)-1, windowRect(3),     windowRect(4)];
    
    HideCursor;
    % raise priority level:
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
    % Query the screen refresh rate:
    ifi = Screen('GetFlipInterval',win); %in sec
    
    % Set the alpha-blending:
    % We want a linear superposition of the dots should they overlap:
    % Just like the Gabors in GarboriumDemo.m (see there for further info).
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    % not sure how this will conflict with the above command
    % about the linear superposition of dots... but it doesn't seem to cause problems
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Make the textures for the dots, 1 for black and white
    Dots(1) = Screen('MakeTexture',win,env,[],[],2);  % white dot
    Dots(2) = Screen('MakeTexture',win,env2,[],[],2); % black dot
    
    % Generate the annulus texture:
    AnnImages = 0.5.*ones(imsize2,imsize2,2); %specify RGB matrix of annulus, grey
    AnnImages(:,:,2) = 1.*J; %specify Alpha channel of annulus
    annulus = Screen('MakeTexture',win,AnnImages,[],[],2);
    
    % Generate the Inner ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaInner;
    fixationRingTextureInner = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Generate the Outer ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaOuter;
    fixationRingTextureOuter = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Preallocate array with destination rectangles:
    % This also defines initial dot locations
    % for the very first drawn stimulus frame:
    texrect = Screen('Rect', Dots(1));
    parameters.TextureRect = texrect; %Save into parameters
    
    %%
    % Now, load some of the pre-genereated stimuli into memory to define the stimulus variables.
    % This is determined by parameters.ConditionOrder, but enter 1 for CD and 3 for IOVD
    % Enter value of zero for P, the probe number, so that stimulus
    % parameters can be saved at setup but not during the experiment.
    % Load CD
    [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR] = SetUpDotPositions ...
        (1, Dots, parameters, 0);
    % Load IOVD
    [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR] = SetUpDotPositions ...
        (3, Dots, parameters, 0);
    
    % Save the parameters just established:
    save(subjectData.filename, 'subjectData', 'parameters');
    
    %%
    %%%%-----------------------------------------%%%%
    %               Welcome screen
    %%%%-----------------------------------------%%%%
    
    % Display the welcome screen and wait for the user to begin.
    Screen('TextFont',win, 'Arial');
    Screen('TextSize',win, 24);
    if useHardwareStereo
        Screen('SelectStereoDrawBuffer', win, 0);   % flag of 0= left eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 8 '' if the 3D motion moves AWAY from you,', ...
            '\n \nor press '' 2 '' if the 3D motion moves TOWARDS you.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        
        Screen('SelectStereoDrawBuffer', win, 1);   % flag of 1= right eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 8 '' if the 3D motion moves AWAY from you,', ...
            '\n \nor press '' 2 '' if the 3D motion moves TOWARDS you.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        Screen('Flip', win); % , [], [], [], 1);
    else
        
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 8 '' if the 3D motion moves AWAY from you,', ...
            '\n \nor press '' 2 '' if the 3D motion moves TOWARDS you.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        Screen('Flip', win); %, [], [], [], 1);
    end
    
    WaitSecs(0.2)
    KbCheck(); % take a quick KbCheck to load it now & flush any stored events
    
    % Wait for user response to continue...
    ButtonPressed = 0;
    while ~ButtonPressed
        % if 'q' is pressed, abort
        [KeyIsDown, ~, keyCode] = KbCheck();
        if KeyIsDown % a key has been pressed
            if keyCode(RespQuit)
                ForcedQuit = true
                ExitGracefully(UsingVP, ForcedQuit)
                %if any other button on the keyboard has been pressed
            else
                ButtonPressed = 1;
            end
        end
    end
    
    % Blank the screen and wait 2 secs before beginning.
    Screen('Flip', win);
    WaitSecs(2);
    missedFrames = 0; % Set up counter for missed frames.
    
    %%%%---------------------------------------------------------------%%%%
    %                     Present the stimuli!
    %%%%---------------------------------------------------------------%%%%
    
    if ~ismac
        jheapcl; %clear the java heap space.
    end
    
    % Set left and right eye contrasts for first few blank frames
    LeftEyeContrast = 0;
    RightEyeContrast = 0;
    
    F = 1; % This value increments across PER EYE frames within a single stimulus (ie. every 2nd video frame) and indexes stimulus variables.
    % It is reset to 1 for each new blank or probe.
    P = 1; % counter for the total probes/blanks
    overallPerEyeFrame = 1; %This increments across all per-eye frames across the experiment, to index the button press matrix. 
    
    if UsingVP
        % ***** send trigger ***
        % we will need a trigger here to mark beginning of stimulus duty cycle...
        Datapixx('SetDoutValues', transformindex(10)); % 10 to indicate start/end of run
        Datapixx('RegWrRd');
    end
    
    % Draw the rings again for 1 frame so that there is no 'blank'
    % frame presented in between, & hopefully stop the 'flickering' of the goggles between intervals
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % need to flip the alpha around again (anti-aliasing)
    Screen('DrawTexture',win, fixationRingTextureInner);
    Screen('DrawTexture',win, fixationRingTextureOuter);
    % Draw the fixation cross:
    Screen('FillRect',win,[0 0 0],fixationCross);
    Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
    Screen('BlendFunction', win, GL_ONE, GL_ZERO);
    Screen('DrawingFinished', win);
    vbl = Screen('Flip',win); %sync vbl to start time
    startTime = vbl;
    
    % For blank periods we just set the contrast to zero and continue presenting the previous set of dots.
    % For each new event we update the dot positions and reset contrast to 1.
    
    while P <= length(AssignCond) % loop through all probes, including blanks
        
        %Determine what the correct length of the stimulus or blank should be, if this is the first frame:
        if F == 1
            probeStartVBL = vbl; %take start point as most recent vbl
            probeEndVBL = probeStartVBL + stimDurations(P) - ifi;
            
            % Send trigger to mark first frame of stimulus:
            if UsingVP
                Datapixx('SetDoutValues', transformindex(AssignCond(P)));
                Datapixx('RegWrRd');
            end
        end
        
        switch AssignCond(P)
            case 5 % it's a blank period
                % Set contrast back to 0 for blank frames.
                LeftEyeContrast = 0;
                RightEyeContrast = 0;
                
            case {1, 2, 3, 4} % It's any of the 4 stimulus categories
                
                %load and set up dot positions
                [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR] = SetUpDotPositions ...
                    (AssignCond(P), Dots, parameters, P);
                
                % Reset dot contrasts back to 1:
                LeftEyeContrast = 1;
                RightEyeContrast = 1;
                % Set aside record of stimulus on its first frame.
                Responses(overallPerEyeFrame,1) = AssignCond(P);
        end
  
        % Draw each stimulus event:
        probeEnd = false;
        while ~probeEnd %This should keep iterating across stim frames until vbl >= probeEndVBL
            
            % Select left-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 0);
            end
            
            %%%%------------------------------------------------%%%%
            %               Draw left eye stimulus:
            %%%%------------------------------------------------%%%%
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('DrawTextures',win, ...
                DotsIdxL(:,mod(F,parameters.FramesFullCycle)+1), ...      % determines the dot texture (colour) to be drawn: ***NOTE: now indexed every frame***
                [], ...                                                   % texture subpart (not using)
                dstRectsL(:,:,mod(F,parameters.FramesFullCycle)+1), ... % determines dot locations for each frame
                [], ...                                                   % rotate texture (not using)
                [], LeftEyeContrast);                                     % Final argument is for contrast (modulates global alpha)
            
            %Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            % Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            if DrawRings
                Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % need to flip the alpha around again (anti-aliasing)
                Screen('DrawTexture',win, fixationRingTextureInner);
                Screen('DrawTexture',win, fixationRingTextureOuter);
                % Draw the black fixation cross:
                Screen('FillRect',win,[0 0 0],fixationCross);
            end
            
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            %Draw blue lines:
            Screen('FillRect', win, [0, 0, 1], blueRectLeftOn);
            Screen('FillRect', win, [0, 0, 0], blueRectLeftOff);
            
            % Select right-eye image buffer for drawing:
            if useHardwareStereo
                Screen('SelectStereoDrawBuffer', win, 1);
            else
                Screen('DrawingFinished', win);
                [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); % update display on next refresh (& provide deadline)
                %overallFrame = overallFrame + 1;
                if missed > 0
                    missedFrames = missedFrames + 1;
                end
            end
            
            %%%%------------------------------------------------%%%%
            %               Draw right eye stimulus:
            %%%%------------------------------------------------%%%%
            Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
            Screen('DrawTextures',win, ...
                DotsIdxR(:,mod(F,parameters.FramesFullCycle)+1), ...     % determines the dot texture (colour) to be drawn: ***NOTE: now indexed every frame***
                [], ...                                                  % texture subpart (not using)
                dstRectsR(:,:,mod(F,parameters.FramesFullCycle)+1), ...  % determines dot locations for each frame
                [], ...                                                  % Rotate texture (not using)
                [], RightEyeContrast);                                   % Final argument is contrast (modulates global alpha)
            
            %Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            % Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            if DrawRings
                Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
                Screen('DrawTexture',win, fixationRingTextureInner);
                Screen('DrawTexture',win, fixationRingTextureOuter);
                %Draw the fixation cross:
                Screen('FillRect',win,[0 0 0],fixationCross);
            end
            
            Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
            Screen('BlendFunction', win, GL_ONE, GL_ZERO);
            %Draw blue lines:
            Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
            Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
            
            Screen('DrawingFinished', win);
            
            [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
            
            %keep record of any missed frames:
            if missed > 0
                missedFrames = missedFrames + 1;
            end
            
            %%%%------------------------------------------------%%%%
            %               Check for button presses
            %%%%------------------------------------------------%%%%
            % Check for valid keyboard response
            [keyIsDown, ~, keyCode] = KbCheck([],scanList);
            if keyIsDown %If they've pressed a button on keypad,
                % determine if it is valid and proceed
                
                % If they press 2 or 8 for a type of motion, store this response at the given per-eye frame.
                % This gives a rough approximation of reaction time.
                % Correctness of the response will need to be inferred offline on the basis of the most recent probe stimulus.
                % Correct = 1, incorrect = 0.
                if any(keyCode(TowardsMotion))            % Respond TOWARDS
                    Responses(overallPerEyeFrame,2) = -1; % -1 to indicate towards response.
                    
                elseif any(keyCode(AwayMotion))           % Respond AWAY
                    Responses(overallPerEyeFrame,2) = 1;  % +1 to indicate away response.
                    
                elseif keyCode(RespQuit)
                    % break out of program if 'q' is pressed
                    ForcedQuit = true
                    ExitGracefully(UsingVP, ForcedQuit)
                end
            end
            
            %%%%------------------------------------------------%%%%
            %               Prepare for next frame:
            %%%%------------------------------------------------%%%%
            % Increment frames. These both increment on a PER-EYE basis
            F = F + 1; % Only F is reset with each probe/blank
            overallPerEyeFrame = overallPerEyeFrame + 1; % to index the button presses
            
            % This is the important bit:
            % If we've reached the determined end time of the probe/blank
            % we reset F and move on to the next one.
            % This means no stimulus ever exceeds its deadline (by more than 1 frame anyway)
            % and we shouldn't miss any either.
            if vbl >= probeEndVBL
                F = 1;           % Reset F, which increments at the PER EYE frame rate
                P = P+1;         % Increment the blanks or probes
                probeEnd = true; % This should terminate the current stimulus execution and move on to next one
            end
            
        end % end of loop across stimulus frames
        
    end  % end of loop across total probes/blanks
    
    % Send final trigger to mark end of stimulus run as '10':
    if UsingVP
        Datapixx('SetDoutValues', transformindex(10));
        Datapixx('RegWrRd');
    end
    
    % Print out some important values.
    missedFrames
    P
    TotalEvents = length(AssignCond)
    runTime = vbl - startTime
    
catch MException
    
    % We throw the error again so the user sees the error description.
    rethrow (MException)
    psychrethrow(psychlasterror);
    ExitGracefully (UsingVP, ForcedQuit)
    error('Error!')
    
end % End of try/catch statement

%%
%%%%-----------------------------------------------------------%%%%
%                       Save the results
%%%%-----------------------------------------------------------%%%%

% Store the data and some info about the run
subjectData.Responses = Responses; % button press data
subjectData.missedFrames = missedFrames;

% Save the data and the updated parameters information:
save(subjectData.filename, 'subjectData', 'parameters')

% Close everything down & exit:
ExitGracefully (UsingVP, ForcedQuit)
end %end of main function

%Now for the sub-functions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExitGracefully (UsingVP, ForcedQuit)
%...need to shut everything down here...

% turn off the prioritisation:
Priority( 0 ); % restore priority

if UsingVP        % close down the ViewPixx or ProPixx
    Datapixx('DisableVideoScanningBacklight');
    if Datapixx('IsViewpixx3D')
        Datapixx('DisableVideoLcd3D60Hz');
    end
    Datapixx('RegWr');
    %Datapixx('Close'); %closing it here might cause it to crash?
end

% Close down the screen:
Screen('CloseAll')

Datapixx('Close'); % closing the Datapixx here (after closing the screen) might stop it from crashing

% Bring back the mouse cursor:
ShowCursor();

% announce to cmd window if the program was aborted by the user
if ForcedQuit
    error('You quit the program!')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR] = SetUpDotPositions (Probe, Dots, parameters, P)

% This is where the pre-generated dot positions are loaded and the disparity added.
% (see GenerateCDdots_Attn_fMRI.m)
% (see GenerateIOVDdots_2.m & GenerateCDdots.m)
% These are picked at random (with replacement) from those already in the folder.
%%
% Indicate how many existing sets of dot stimuli are in the subfolder 'stimuli' for loading later.
StimInFolder = 5;

% Convert frequency into a string so we can load the appropriate stimulus:
IOVDstrfreq = num2str(parameters.IOVD_freq);
IOVDstrfreq(IOVDstrfreq == '.') = '_'; % replace decimal point with underscore
CDstrfreq = num2str(parameters.CD_freq);
CDstrfreq(CDstrfreq == '.') = '_';

% Draw a random value to load the stimuli with:
RandDraw = ceil(rand*StimInFolder);

% Set up binary control for stimuli.
% This single matrix will determine the type of motion and the direction
% Column 1: 0 for CD, 1 for IOVD.
% Column 2: 1 for towards, 0 for away.
MotionSwitch = [0 1; 0 0; 1 1; 1 0];

% Determine the direction of the motion.
% The value of -1 or 1 is multiplied by the sawtooth wave of the motion to give us
% towards or away motion.
if MotionSwitch(Probe, 2) % value of 1 == towards motion.
    DirnConstant = -1; % a negative sawtooth wave gives us TOWARDS motion
else
    DirnConstant = 1; % a positive sawtooth wave gives us AWAY motion
end

%%
% *** ---------------------------------- ***
%               Set up CD
% *** ---------------------------------- ***
if MotionSwitch(Probe, 1) == 0
    
    % Load random set of dot positions
    FileNameCD = fullfile('stimuli',['CD_dots_', CDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
    CD_stim = load(FileNameCD); % includes parameters file and dot positions matrices
    % Just cut the parameters back to a smaller struct:
    CDparams = CD_stim.DotParams;
    
    % Assign the pre-generated dot position matrices.
    % First, we duplicate the left eye dot positions to make the right eye dot positions.
    % At this stage all dots are identical for the two eyes (disparity has not been added yet).
    dot_posL = CD_stim.dot_posL;
    dot_posR = CD_stim.dot_posL;
    
    % Set aside the number of frames in a single cycle in parameters,
    % important for later stimulus control
    parameters.FramesFullCycle = CDparams.FramesFullCycle;
    
    % set aside the disparity
    CDparams.disparity = parameters.CD_disparityPixels;
    
    % Set aside the frequency (this is not stored when stimuli are generated)
    CDparams.frequency = parameters.CD_freq; %(defined above)
    
    % set up a few extra parameters of the sine waves for each cue:
    % the period, in sec
    CDparams.period = 1/CDparams.frequency;
    
    % the angular frequency of the sine waves:
    CDparams.angFreq = 2 * pi * CDparams.frequency;
    
    % Length of the sine wave:
    % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
    % end of the cycle and the beginning of the next cycle.
    % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
    %CDparams.t = linspace(CDparams.period/CDparams.FramesFullCycle, CDparams.period, CDparams.FramesFullCycle);
    % *** NOTE: *** not so with the sawtooth wave. We don't want the wave to wrap at this stage, so just one steady change in disparity.
    % So we start at zero and go one frame less than the period.
    CDparams.t = linspace(0, CDparams.period - (CDparams.period/CDparams.FramesFullCycle), CDparams.FramesFullCycle);
    
    
    % Now make one full cycle of the SAWTOOTH wave, no matter the frequency
    % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
    % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
    % *** Here, we also multiply by the DirnConstant to determine whether the motion is towards or away (as defined above).
    CDparams.SineWave =  DirnConstant * CDparams.disparity * sawtooth(CDparams.angFreq * CDparams.t); %disparity * sin(angFreq * t)
    
    % assign dstRects: destination rect matrices
    dstRectsL = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
    dstRectsR = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
    
    % Shift dot trajectory: add disparity.
    % We do this for every frame in the experiment, in units of integer (per-eye) frames!
    for fr = 1:CDparams.FramesFullCycle
        % determine dot coordinates: remember, y position does not change: horiz disparity only
        % Update dot position according to sawtooth trajectory on each (per-eye) frame
        % Left eye: -sin
        dstRectsL(:,:,fr) = CenterRectOnPoint(parameters.TextureRect, ...   % the size of the dot texture
            dot_posL(:,1,fr) - CDparams.SineWave(fr), ...                   % the x positions
            dot_posL(:,2,fr))';                                             % the y positions
        % Right eye: +sin
        dstRectsR(:,:,fr) = CenterRectOnPoint(parameters.TextureRect, ...
            dot_posR(:,1,fr) + CDparams.SineWave(fr), ...
            dot_posR(:,2,fr))';
    end
    
    % Now we have the dot position indices (the dstRects), define the dot texture indices.
    % These are simply half white (Dots(1) and half black (Dots(2)) in each eye.
    
    % NOTE: since fixing the Right eye IOVD flicker bug (March 2016), DotsIdx matrices are now 2D with num_dots * FramesFullCycle dimensions.
    % So we need to set up the matrix accordingly. We won't do it like the below any more:
    % DotsIdxL = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
    % DotsIdxR = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
    
    % Because in CD the dots only last one frame & then move, it doesn't really matter how their colour is assigned (as long as there is equivalent
    % numbers of black and white AND they are the same across the eyes.
    % So we can set them up in the same way as IOVD, except they are not later sorted by colour (using MinDist), and they are identical in the 2 eyes:
    
    %Left eye:
    DotsIdxL = [repmat(Dots(1),1,CDparams.NumDots/2); repmat(Dots(2),1,CDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
    DotsIdxL = repmat(reshape(DotsIdxL,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle);   % Replicate for each frame
    %Right eye:
    DotsIdxR = [repmat(Dots(1),1,CDparams.NumDots/2); repmat(Dots(2),1,CDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
    DotsIdxR = repmat(reshape(DotsIdxR,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle);   % Replicate for each frame
    
    % Set aside details of the stimulus parameters for this PROBE.
    % Most of these are already stored in the pre-generated stimulus files,
    % but it is probably worthwhile to have them here again in a single output file
    
    % NOTE: storing the file name of each probe was potentially causing
    % missed frames, so we will just have to do without it.
    if P == 0
        parameters.CD_parameters = CDparams; % Store the CD parameters once only, during set-up
        %     else
        %         parameters.RandomStimulusNoEachProbe{P} = FileNameCD;
    end
    
    % *** ---------------------------------- ***
    %               Set up IOVD
    % *** ---------------------------------- ***
    
else % if MotionSwitch(Probe, 1) ~= 0, then it is IOVD.
    
    FileNameIOVD = fullfile('stimuli',['IOVD_dots_', IOVDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
    IOVD_stim = load(FileNameIOVD); % includes parameters file and dot positions matrices
    % Just cut the parameters back to a smaller struct:
    IOVDparams = IOVD_stim.DotParams;
    
    % Set aside the number of frames in a single cycle in parameters,
    % important for later stimulus control
    parameters.FramesFullCycle = IOVDparams.FramesFullCycle;
    
    % set aside the disparity
    IOVDparams.disparity = parameters.IOVD_disparityPixels;
    % For reference also store the disparity in Arcmin
    IOVDparams.DisparityInArcmin = (parameters.IOVD_disparityPixels/parameters.PPD)*60;
    
    % Set aside the frequency (this is not stored when stimuli are generated)
    IOVDparams.frequency = parameters.IOVD_freq; %(defined above)
    
    % set up a few extra parameters of the sine waves for each cue:
    % the period, in sec
    IOVDparams.period = 1/IOVDparams.frequency;
    
    % the angular frequency of the sine waves:
    IOVDparams.angFreq = 2 * pi * IOVDparams.frequency;
    
    % Length of the SAWTOOTH wave:
    % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
    % end of the cycle and the beginning of the next cycle.
    % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
    IOVDparams.t = linspace(IOVDparams.period/IOVDparams.FramesFullCycle, IOVDparams.period, IOVDparams.FramesFullCycle);
    
    % Now make one full cycle of the SAWTOOTH wave, no matter the frequency
    % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
    % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
    % *** Here, we also multiply by the DirnConstant to determine whether the motion is towards or away (as defined above).
    IOVDparams.SineWave =  DirnConstant * IOVDparams.disparity * sawtooth(IOVDparams.angFreq * IOVDparams.t); %disparity * sin(angFreq * t)
    
    % assign dstRects: destination rect matrices
    dstRectsL = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
    dstRectsR = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
    
    % Shift dot trajectory: add disparity.
    % Remember that with IOVD the dot positions for each eye are NOT identical (before the horizontal shifts are added)
    % as they are with CD, that is why we have 2 matrices, 1 for each of the 2 eyes
    % Dot lifetime is determined by how many dot positions are re-randomised with each frame (see GenerateIOVDdots.m)
    
    % Determine a value to modulate the sine wave by.
    % In this way the sinusoidal disparity is either subtracted or added to the dot positions with each frame.
    % If both eyes are added (as in the IOVD_control), then the dots move in the same direction in the 2 eyes.
    % The two columns in SineModulate are for all dots in each eye: col 1 for LE, col 2 for RE.
    
    % *** In this paradigm (MID_SSVEP) we probably won't be using the control stimulus so there is only 1 value for 'SineModulate'.
    % It's probably therefore not necessary but we will leave it in in case later we
    % turn around and decide to display the IOVD control stimulus.
    
    % subtract left eye, add right eye for MID.
    % We will assign this variable here, but if we are using the IOVD control as our 'probe'
    % (rather than simple noise), we will need to re-assign it anyway
    SineModulate = ones(IOVDparams.NumDots,2);
    SineModulate(:,1) = -1; %-1 * sin for the left eye only in actual IOVD; +1 * sin for right eye
    
    % Set up another SineModulate vector, of the type that would be used if we were running the IOVD control.
    % Here, 50% of dots will go left (of both colours), 50% will go right
    %SineModulate2 =  ones(IOVDparams.NumDots,2);
    %SineModulate2(1:end/2,:) = -1;
    
    for fr = 1:IOVDparams.FramesFullCycle
        % determine dot coordinates: remember, y position does not change: horiz disparity only
        % Update dot position according to sinsoidal trajectory on each (per-eye) frame
        % Left eye: -sin
        dstRectsL(:,:,fr) = CenterRectOnPoint(parameters.TextureRect, ...   % the size of the dot texture
            IOVD_stim.DotPosBinoc{1}(:,1,fr) + ...                          % original positions
            SineModulate(:,1) * ...                                         % sign of the sine wave, determines direction of shift (important if we want to use the IOVD control at some point).
            IOVDparams.SineWave(fr), ...                                    % add the shift to the x positions
            IOVD_stim.DotPosBinoc{1}(:,2,fr))';                             % y positions do not change.
        
        %Right eye: +sin
        dstRectsR(:,:,fr) = CenterRectOnPoint(parameters.TextureRect, ...
            IOVD_stim.DotPosBinoc{2}(:,1,fr) + ...
            SineModulate(:,2) * ...
            IOVDparams.SineWave(fr), ...
            IOVD_stim.DotPosBinoc{2}(:,2,fr))';
    end
    
    % Now we have the dot position indices (the dstRects), define the dot texture indices.
    % These are simply half white (Dots(1) and half black (Dots(2)) in each eye, same as CD.
    
    % *** This is where the dot color must be assigned in an alternate fashion for L & R eye
    % to make (mostly) sure that neighbouring dots across the eyes are of a different color
    % However you do it, colour must be assigned oppositely to the 2 eyes for this to work
    
    % Previously we only assigned dot color in alternating fashion for IOVD, not the control.
    % This was because we previously used the same set of dots in the 2 eyes for the control,
    % so if color was assigned in an alternating fashion then the dots would flicker between the 2 colours as the
    % dots were presented in alternating fashion to the 2 eyes.
    % Now, because we use a different set of dots, it is ok for the dots to be assigned in alternating manner
    % across the eyes. This also maintains the method used to reduce binocular correlations in the new IOVD algorithm, applying it to the control also.
    
    % *** To solve the R eye flicker bug associated with colours changing across dot lifetime, the DotIdx matrices
    % must now be a 2D matrix with FrmsFullCycle * num_dots dimensions
    % Hence, both DotsIdxL & DotsIdxR must now be indexed by the current frame in the cycle as the dots are drawn
    % (as was done with the dot position matrices). ***
    
    % Generate the alternating matrices of dot colours. Note that they need to alternate.
    % If we just assigned top half of the matrix to black, bottom half white (for eg.), then that would mean
    % mostly white dots would go one direction in the IOVD control stimulus, mostly black in the other direction, for example.
    %(this is because dots are now positioned in alternating strips across the 2 eyes)
    % We can't have that; in the control we would want equal numbers of white/black dots going in either direction.
    
    % Left eye: WBWBWB ....
    DotsIdxL = [repmat(Dots(1),1,IOVDparams.NumDots/2); repmat(Dots(2),1,IOVDparams.NumDots/2)]; % Dots(1) = white; Dots(2)=black
    DotsIdxL = repmat(reshape(DotsIdxL,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % *Replicate for every per-eye frame in the whole run*
    % Right eye: BWBWBW....
    DotsIdxR = [repmat(Dots(2),1,IOVDparams.NumDots/2); repmat(Dots(1),1,IOVDparams.NumDots/2)];
    DotsIdxR = repmat(reshape(DotsIdxR,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % *Replicate for every per-eye frame in the whole run*
    
    % Now sort these dot colour matrices according to the indices established earlier that account for dot lifetime
    % Because dot colour was assigned in an opposite manner above, we sort the matrices using the same indices.
    % This retains the opposite polarity assignment across eyes
    DotsIdxR = DotsIdxR(IOVD_stim.MinDist);
    DotsIdxL = DotsIdxL(IOVD_stim.MinDist);
    
    % Set aside details of the stimulus parameters for this PROBE.
    % Most of these are already stored in the pre-generated stimulus files,
    % but it is probably worthwhile to have them here again in a single output file
    
    % NOTE: storing the file name of each probe was potentially causing
    % missed frames, so we will just have to do without it.
    if P == 0
        parameters.IOVD_parameters = IOVDparams; % Store the IOVD parameters once, during set-up
        %     else
        %         parameters.RandomStimulusNoEachProbe{P} = FileNameIOVD;
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rect = fixation_cross(width,height,centrex,centrey)
% $Id: fixation_cross.m 4 2007-06-23 11:13:19Z damienm $
%Make a small fixation cross to display on screen.
%Width and height are in pixels.
%centrex & centrey give the x & y screen coordinates where you want the cross centred.

rect = zeros(4,2);

width = width/2;
height = height/2;

rect(1,1) = -height;
rect(2,1) = width;
rect(3,1) = height;
rect(4,1) = -width;

rect(1,2) = -width;
rect(2,2) = height;
rect(3,2) = width;
rect(4,2) = -height;


rect(1:2:4,:) = rect(1:2:4,:) + centrex;
rect(2:2:4,:) = rect(2:2:4,:) + centrey;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = transformindex(input)

% fixes the binary inputs for the EEG amplifier because the pins are in a different order from the ViewPixx
% desired numbers must be <256
% DHB 18/8/14

truebits = 2.^(2:2:24);
dn = dec2bin(input,length(truebits));
output = 0;
for m = 1:length(truebits)
    output = output + truebits(m)*str2num(dn(end-m+1));
end

end