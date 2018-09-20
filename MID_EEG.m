function MID_EEG (SubjCode, Cue, direction, runN)

% Frequency should be in Hz
% Amplitude should be in ARCMIN!
% These can be defined below.
% direction = 'towards' or 'away', as a string. 

% towards: negative sawtooth wave
% away: positive sawtooth wave.

% **** sawtooth wave should be ok! ***

%%
% Define some initial parameters.

%%%%---------------------------------------%%%%
% Specific for S-cone and achromatic thresholds:

% Colour parameters:
UseS_cone = true; %Flag for S cone stimulus or L-M cone stimulus (=1 for S cone)

% If you want the 2 deg cone fundamentals, set to true:
Load2Deg = false; % otherwise, the 10 deg ones will be loaded.

% Testing S-cone thresholds or black and white?
TestScone = false; % true for s-cone, false for black/white thresholds.

% Contrast of the achromatic dots.
% Apparently this should be only about 10% of the cone contrast, but specify here:
AchromContrast = 1.0; 

%%%%---------------------------------------%%%%
% Others:

% Define some of the display parameters:
PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)

% If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = false;

% Switch for whether you want the annulus superimposed over the dots:
DrawAnnulus = true;

% Draw the fixation lock/rings
DrawRings = true;

% Spatiotemporal parameters for IOVD and CD:
IOVD_freq = 4;
CD_freq = 4;

% Define the max. change in disparity (horiz shift of dots) for the 2 cues:
% the numerator is in arcmin (but converted to pixels) akin to the amplitude of the sine wave
IOVD_disparity = 128/60 * PPD;
CD_disparity = 32/60 * PPD;

%%
%%%%---------------------------------------%%%%
%       create a unique file name
%%%%---------------------------------------%%%%

% Colour (s-cone) or achromatic thresholds?
if TestScone
    ColourType = 'Scone';
else
    ColourType = 'Achrom';
end

% Set up file name
% Enter the frequency and amplitude into file name
fName = fullfile('data', [SubjCode, '_', 'MID_EEG_', Cue, '_', direction, '_', ColourType, '_', ...
    num2str(runN), '.mat']);

%Check that it doesn't already exist
if exist( fName, 'file' )
    userResponse = input( 'WARNING: File exists. Overwrite? Enter y or n: ', 's' );
    if ~strcmp( userResponse, 'y' )
        error('Aborting function' );
    end
end

%%
%Choose the screen: it is usually the max screen no. available.
%Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
%So if we're on that machine, need to -1 from the screen number:
[~, CompName] = system('hostname'); %find out the computer name
if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ... %normal strcmp not working here, can't figure out why...
        && (length(Screen('Screens'))>1) %and there is more than 1 display connected...
    WhichScreen = max( Screen( 'Screens' ) )-1;
else
    WhichScreen = max( Screen( 'Screens' ) ); %should be right for any other machine!
end

screenRect = Screen('Rect',WhichScreen); %get the screen resolution.
centreX = screenRect(3)/2;
centreY = screenRect(4)/2;
RefreshRate = Screen('NominalFrameRate', WhichScreen);

% Determine the file name that contains the theta values for isoluminance.
% These must have been specified previously.
% NOTE: Here we assume that you are using the first run of the isolum task (_1).
% It would be unusual if you were using another run...
% AND we are assuming it's for S cone. We can build in other functionality there later if needed.
subjectData.IsoFileName = fullfile('data', ...
    [SubjCode, '_isolum_S_cone_1.mat']);

%%%%-------------------------------%%%%
%         Load the spectra:           %
%%%%-------------------------------%%%%

% If we're not on a Vpixx device (or not using it), we just load
% the Viewpixx spectra anyway.
% This task should only be performed on the Viewpixx 3D.
% It is easier to load these here so they can be parsed to the LMS2RGB function.
% We also load the subject's isoluminance values here too.
if UsingVP
    Datapixx('Open');
    if Datapixx('IsViewpixx3D')  % If it's the Viewpixx3D
        load('C:\Users\PS-RMLab\Documents\Calibration\Viewpixx_Processed_cal_data_25_8_2016.mat')
        % Load the isoluminance values:
        load(subjectData.IsoFileName);
    end
else % If we're not using a Vpixx device, just load the Viewpixx spectra anyway (need spectra from somewhere)...
    load('C:\Users\PS-RMLab\Documents\Calibration\Viewpixx_Processed_cal_data_25_8_2016.mat')
    % Load the isoluminance values:
    load(subjectData.IsoFileName);
end

% Save the information about the PC hardware:
subjectData.HostPC = CompName;

if ~ismac
    jheapcl; %clear the java heap space.
end

%%
%%%%-------------------------------%%%%
%       define response keys
%%%%-------------------------------%%%%

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

% Subjects can press '1' for signal in interval 1: either on the keyboard or the 1 on the num keypad are fine.
Interval1 = [KbName('1!'), KbName('1')];

% Or press '2' for signal in interval 2.
Interval2 = [KbName('2@'), KbName('2')];
% Finally, press a q at any time to quit the program.
RespQuit = KbName('q'); % q to quit.

% Add the keycodes for the appropriate buttons on the gamepad, if we're
% using the mac setup.
% Now, at time of writing we are using the 'Afterglow gamepad for PS3' device.
% We want to use the L2 and R2 buttons on the top left & right of the pad,
% for the 1st and 2nd intervals, respectively.
% Because the gamepad registers as a HID device, like a keyboard, we can
% define the button responses in the same way as we would a keyboard response.
% On the Afterglow gamepad, L2 = button 7; R2 = button 8.
% They are equivalent to KbName('d') & KbName('e'), on the keyboard, respectively.
% This of course means that if a subject presses 'd' or 'e' on the keyboard
% (as well as 1 or 2), that will register as a response, so we will just
% have to warn against that.

if ismac
    Interval1 = [Interval1, KbName('d')];
    Interval2 = [Interval2, KbName('e')];
    % Also get the Gamepad indices:
    GamepadIndex = GetGamepadIndices;
end

%%
%%%%-------------------------------%%%%
%           Set up stimulus
%%%%-------------------------------%%%%

% Dot colour parameters:

% Load the Stockman/Sharpe cone fundamentals:
if Load2Deg % for 2 deg fundamentals
    load('C:\Users\PS-RMLab\Documents\Colour\StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
    %load('StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
else        % or for the 10 deg fundamentals
    load('C:\Users\PS-RMLab\Documents\Colour\StockmanSharpe_10deg_cone_fundamentals_1nm.mat')
    %load('StockmanSharpe_10deg_cone_fundamentals_1nm.mat')
end

% Timing:
TrialDuration = 10; % duration of the stimulus, IN SEC

% the effective frame rate PER EYE: since each eye is stimulated on successive video frames
% Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
PeyeFR = RefreshRate/2; % Per eye f.r. is always half the absolute f.r.

% stimuli for each eye are indexed with the same number: they increment on a per-eye frame rate basis
% But PTB refresh 'flips' of course happen at the machine frame rate: there are 2 flips within a single stimulus loop: one for each eye.
% PeyeFrames is the number of (per-eye) frames needed to achieve the trial duration.
% This is an important number and will determine how many points in the sine wave are sampled (& whether the points are wrapped around).
% Note that 'PeyeFrames' is determined by the stimulus duration and may not be the same as 'FramesFullCycle'
PerEyeFrames = round(PeyeFR*TrialDuration); % also round to integer number of frames

% Contrast:
pdc = 1; % the peak dot contrast

% Define the dot texture, a square-shaped sheet of dots.
% Make the texture the same size as the height of the screen
% (well really we are using many little dot textures, but 'imsize' defines the region they inhabit)
% Note, this should really be the same for both CD and IOVD...
imsize = screenRect(4);

% specify dot size:
dot_sigma_in_degrees = 0.05;            % 0.064; %size of SD of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; % sigma in pixels
dotsize = round(dot_sigma * 10);        % make the dots some multiple of sigma
% NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
% It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.

% RM 12/5/16
% Make the dots a cosine profile, rather than Gaussian.
% This is because the S-cone dots are too hard to see.
circ_dot_size = dot_sigma * 2.355 * 2; % this approximates the full dot width of the gaussian dots, since there are 2.355 sigma under the curve at FWHM
dot_cos_smooth = circ_dot_size / 3; % make the smoothing 1/3 width of the dots

%Set up the cosine edges for the dots:
dot_rad = circ_dot_size/2;
x0 = (dotsize+1)/2;
y0 = (dotsize+1)/2;
env = ones(dotsize);
for (ii=1:dotsize)
    for (jj=1:dotsize)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > dot_rad^2)
            env(ii,jj) = 0;
        elseif (r2 > (dot_rad - dot_cos_smooth)^2)
            env(ii,jj) = cos(pi.*(sqrt(r2)-dot_rad+dot_cos_smooth)/(2*dot_cos_smooth))^2;
        end
    end
end

% % This gives us a white-peaked dot (+ve contrast polarity)
% env = exp(-r.^2/(2*dot_sigma.^2)); % gaussian window
% env = env./max(max(abs(env)));     % normalize peak to +/- 1
img = repmat(env,1,1,3); %pre-allocate the image matrix.
img2 = img; % for the negative polarity dots.

% set up the raised cosine annular window.
% specify parameters for the annulus.
% Note that these are the same as with the fMRI s-cone experiment.
inrad = PPD * 0.5;     % inner radius of annulus (in pixels), for fixation spot
outrad = PPD * 7.5; % outer radius of annulus (in pixels)
% define extent of spatial raised cosine at edge of aperture (in pixels)
cos_smooth = dotsize; % make it one dot size wide
% This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2;   % double the texture size
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

% Set aside some more of the stimulus parameters:
parameters.DotColourType = ColourType;
parameters.Cue = Cue;
parameters.MotionDirection = direction;
parameters.AchromaticContrast = AchromContrast;
parameters.TrialDurationInSec = TrialDuration;
parameters.PerEyeFrames = PerEyeFrames;
parameters.PeakDotContrast = pdc;
parameters.DotSigmaInDeg = dot_sigma_in_degrees;
parameters.OuterStimRadInDeg = outrad/PPD;
parameters.InnerStimRadInDeg = inrad/PPD;

%%
%%%%-------------------------------%%%%
%       Set up fixation
%%%%-------------------------------%%%%

% Set up the fixation cross or spot:
% This is drawn directly to the screen using Screen('FillRect')
crossWidth = 2;
crossHeight = 10;
fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

% Make the fixation lock rings:
% We have an inner one around fixation and an outer one right on the edge of screen.
% These could probably be defined as a single texture (rather than 2) but I think that will complicate matters with the alpha-blending settings.
% (they are complicated enough already)
ringRadiusInner = PPD*0.4;                % ring surrounding fixation
ringRadiusOuter = screenRect(4)/2;        % outer edge (radius) of the ring: the edge of the screen
ringWidthInner = ringRadiusInner - PPD/5; % 1/5 of a degree thick
ringWidthOuter = ringRadiusOuter - PPD/3; % 1/3 of a degree thick

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
ForcedQuit = false; % this is a flag for the exit function to indicate whether the program was aborted

try % Start a try/catch statement, in case something goes awry with the PTB functions
    
    %----------------------------
    %     Set up the screen
    %----------------------------
    
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
    if UsingVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        % The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        Datapixx('EnableVideoStereoBlueline');
        Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses
        % Prepare for the M16 mode for 16 bit luminance output precision: important when varying contrast
        % PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay');
        
        if Datapixx('IsViewpixx3D') % If it's the Viewpixx3D
            
            % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
            % SHOULD NOT use Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
            % The values are the RGB channel gamma function fits as measured through the goggles, 25/8/16.
            % The spectral measurements were made using the Jaz spectrometer on the ShuttlePC
            
            %NOTE: gamma correction might make the dots (esp the black/white dots) seem 'strange',
            % but it is done correctly and the same way as it is done before the isoluminance task so it must be that way (RM: 1/Aug/17).
            
            gammaRed = mean([GammaValues(1,1,1), GammaValues(1,1,2)]);
            gammaGreen = mean([GammaValues(2,1,1), GammaValues(2,1,2)]);
            gammaBlue = mean([GammaValues(3,1,1), GammaValues(3,1,2)]);
            
            % We'll use the average of the right and left gamma values
            PsychColorCorrection('SetEncodingGamma', win, [1/gammaRed, 1/gammaGreen, 1/gammaBlue]);
            
            % Just use the grey measurements
            %PsychColorCorrection('SetEncodingGamma', win, 1/mean([GammaValues(4,1,1), GammaValues(4,1,2)]));
            
            %Datapixx('EnableVideoLcd3D60Hz');
            Datapixx('DisableVideoLcd3D60Hz'); %=> weirdly, this seems to cause less crosstalk (see DBaker)
            subjectData.DisplayType = 'Viewpixx3D'; % set aside the device type for reference
            Datapixx('RegWr');
            
        elseif Datapixx('IsPropixx') % if it's the Propixx DLP projector
            
            subjectData.DisplayType = 'PROpixx';        % set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram',0); % set to normal RGB video processing for driving the LEDs & DLP MMDs
            %Datapixx('RegWr');
            Datapixx('RegWrRd');                        % seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        end
    end
    % Define the 'blue line' parameters
    blueRectLeftOn   = [0,                 windowRect(4)-1, windowRect(3)/4,   windowRect(4)];
    blueRectLeftOff  = [windowRect(3)/4,   windowRect(4)-1, windowRect(3),     windowRect(4)];
    blueRectRightOn  = [0,                 windowRect(4)-1, windowRect(3)*3/4, windowRect(4)];
    blueRectRightOff = [windowRect(3)*3/4, windowRect(4)-1, windowRect(3),     windowRect(4)];
    % No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
    HideCursor;
    % raise priority level:
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
    % Query the screen refresh rate:
    ifi = Screen('GetFlipInterval',win); % in sec
    
    % Set the alpha-blending:
    % We want a linear superposition of the dots should they overlap:
    % Just like the Gabors in GarboriumDemo.m (see there for further info).
    % *** this should be switched off here for the colour to work ***
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    % not sure how this will conflict with the above command
    % about the linear superposition of dots... but it doesn't seem to cause problems
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    %%
    % Make the textures for the dots.
    % Work out the RGB values for each type of dot.
    % We have 6 kinds of dots: black and white,
    % violet/yellow left eye and violet/yellow right eye.
    % To make life easier we will just run through this 6 times
    
    thetas = Col_results.MeanThetas; % Mean thetas for left and right eyes.
    LMS = [1 1 1]; % For the white dot (both eyes)
    
    % Place LMS values in for left and right eyes
    if UseS_cone
        LMS(2,:) = [cos(thetas(1))/sqrt(2), cos(thetas(1))/sqrt(2), sin(thetas(1))]; %LE
        LMS(3,:) = [cos(thetas(2))/sqrt(2), cos(thetas(2))/sqrt(2), sin(thetas(2))]; %RE
    else
        LMS(2,:) = [-sin(thetas(1)), cos(thetas(1)), 0];
        LMS(3,:) = [-sin(thetas(2)), cos(thetas(2)), 0];
    end
    
    % Go through each pixel in the dot and assign appropriate RGB values using LMS2RGB_Vpixx
    % Compute RGB values for every pixel in LMS cone activation space.
    % The first entry below is for the black/white dots' contrast.
    % Although we request full contrast here for the S-cones, we only end up using the max anyway.
    DesiredScale = [ColParams.MaxPossibleContrast * AchromContrast, ColParams.cone_contrast, ColParams.cone_contrast]; %the scale factor (akin to desired cone contrast, I think);
    LMSidx = 0;
    for MakeDots = 1:2:6 % loop across B/W dots, Y/V for L, Y/V for R
        
        % Re-assign the image matrices.
        img = repmat(env,1,1,3);
        img2 = img;
        LMSidx = LMSidx+1; % Counter to index the LMS values for each colour
        stimLMS.dir = LMS(LMSidx,:);
        for ii = 1:dotsize
            for jj = 1:dotsize
                
                %Do positive polarity dot first:
                stimLMS.scale = env(ii,jj) * DesiredScale(LMSidx); %Work out cone activation for this pixel
                [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
                img(ii, jj, :) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
                
                %Now negative polarity dot:
                stimLMS.scale = env(ii,jj) * -DesiredScale(LMSidx); %Note, negative scale for this dot
                [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
                img2(ii, jj, :) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
                
            end
        end
        
        % If you want to make a picture of how the dots look:
%         figure, imshow(img);
%         %put the background in.
%         %Note that currently, we simply use the default background of 0.5, 0.5, 0.5
%         background = [0.5 0.5 0.5];
%         set(gcf, 'color', background);
%         figure, imshow(img2);
%         background = [0.5 0.5 0.5];
%         set(gcf, 'color', background);
        
        %Also plop an alpha channel in to the 4th value of the 3rd dimension.
        %Note that these are identical for the 2 polarity dots and describe the transparency values of each
        %1=opaque, 0=fully transparent.
        img(:,:,4) = env;
        img2(:,:,4) = env;
        
        Dots(MakeDots) = Screen('MakeTexture',win,img,[],[],2);
        Dots(MakeDots+1) = Screen('MakeTexture',win,img2,[],[],2);
    end
    
    % Make the textures for the dots, 1 for black and white
    %Dots(1) = Screen('MakeTexture',win,env,[],[],2);  % white dot
    %Dots(2) = Screen('MakeTexture',win,env2,[],[],2); % black dot
    
    % Generate the annulus texture:
    AnnImages = 0.5.*ones(imsize2,imsize2,2); % specify RGB matrix of annulus, grey
    AnnImages(:,:,2) = 1.*J;                  % specify Alpha channel of annulus
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
    
    % If motion coherence < 1.0, we will need to adjust the dot x positions of the newly-determined noise dots so they are in the centre of the screen.
    % If we don't do this they will be shifted off to the left of screen. This is because imsize is smaller than the full width (ie x distance) of the screen
    AdjustXBy = (screenRect(3) - screenRect(4))/2; % shift X dot positions by this amount to centre them, since image size <screenRect(3)
    
    % Pick and load a stimulus first-off, just to set aside memory for the variables but also to
    % store the cue stimulus parameters now (rather than on each trial).
    % trial_order must be entered as a scalar. Here we enter -1 as a flag to store the parameters & indicate it's not a trial yet
    C = 1;
    Trial = 0;
    [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
        (Cue, -1, C, Trial, amplitude, frequency, AdjustXBy, PPD, Dots, texrect, imsize, parameters, TestScone);
    
    %%
    % Display the welcome screen and wait for the user to begin.
    Screen('TextFont',win, 'Arial');
    Screen('TextSize',win, 24);
    if useHardwareStereo
        Screen('SelectStereoDrawBuffer', win, 0);   % flag of 0= left eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 1 '' if the smooth 3D motion is in the 1st presentation interval,', ...
            '\n \nor press '' 2 '' if the smooth 3D motion is in the 2nd presentation interval.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        
        Screen('SelectStereoDrawBuffer', win, 1);   % flag of 1= right eye
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 1 '' if the smooth 3D motion is in the 1st presentation interval,', ...
            '\n \nor press '' 2 '' if the smooth 3D motion is in the 2nd presentation interval.', ...
            '\n \nPress '' q '' on the keyboard to quit at any time.', ...
            '\n \nPress any button/key to begin.'], ...
            'center', 'center', 0);
        Screen('Flip', win); % , [], [], [], 1);
    else
        
        DrawFormattedText(win,['Welcome ' SubjCode , ...
            '. \n \nPress '' 1 '' if the smooth 3D motion is in the 1st presentation interval,', ...
            '\n \nor press '' 2 '' if the smooth 3D motion is in the 2nd presentation interval.', ...
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
        if ismac
            %  Check both the keyboard and the gamepad;
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
            else %if not, check the response pad
                [ButtonPressed, ~, ~] = KbCheck(GamepadIndex);
            end
        else % if we're not on the mac:
            [KeyIsDown, ~, keyCode] = KbCheck();
            if KeyIsDown
                if keyCode(RespQuit)
                    ForcedQuit = true
                    ExitGracefully(UsingVP, ForcedQuit)
                    %if any other button on the keyboard has been pressed
                else
                    ButtonPressed = 1;
                end
            end
        end
    end
    
    % Blank the screen and wait 1 secs before beginning.
    Screen('Flip', win);
    WaitSecs(1);
    missedFrames = 0;
    
    %%
    %%%%-----------------------------------------------------------%%%%
    %               Perform tasks across trials
    %%%%-----------------------------------------------------------%%%%
    
    for Trial = 1:num_stim  % loop across trials
        
        %---------------------------------------------------------%
        %  Determine current stim based on current Psi staircase
        %---------------------------------------------------------%
        
        % Select which psi to present (random order determined above)...
        trident = psi_order(Trial)
        Trial                       % print out trial number
        stair{trident}.alphaHat     % print out current threshold estimate (NaN on first trial)
        
        % Determine the signal coherence/contrast for the next trial (and print out)
        C = stair{trident}.xLevels( stair{trident}.currentXIndex );
        
        % Here previously we would have set the contrast if measuring contrast thresholds.
        % Since we're measuring coherence thesholds, contrast is at the fixed value set when defining the RGB values for colour or achromatic stimuli.
        TrialContrast = cont;   
        
        % Draw the fixation rings/cross for 1 sec to indicate the start of a trial and prepare the subject:
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % need to flip the alpha around again (anti-aliasing)
        Screen('DrawTexture',win, fixationRingTextureInner);
        Screen('DrawTexture',win, fixationRingTextureOuter);
        % Draw the fixation cross:
        Screen('FillRect',win,[0 0 0],fixationCross);
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); % Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        Screen('DrawingFinished', win);
        Screen('Flip', win);
        WaitSecs(1);
        
        for Interval = 1:2            % Loop across the 2 intervals.
            
            %---------------------------------------------------------%
            %           Load pre-generated dot positions
            %---------------------------------------------------------%
            
            % Determine and load the sequence parameters for the current interval of the current trial.
            % Note that trial_order is entered as a scalar to define the current interval where 0=noise, 1=signal
            % SetUpDotPositions is a nested sub-function defined below
            tic
            %C =1; % if you just want to present a 100% coherent stimulus.
            [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
                (Cue, trial_order(Trial, Interval), C, Trial, amplitude, frequency, AdjustXBy, PPD, Dots, texrect, imsize, parameters, TestScone);
            toc
            
            %%%%---------------------------------------------------------------%%%%
            %                     Present the stimuli!
            %%%%---------------------------------------------------------------%%%%
            
            if ~ismac
                jheapcl; %clear the java heap space.
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
            
            f0 = 1; % This value increments across PER EYE frames (ie. every 2nd video frame): for contrast
            % Also select a random start point within the cycle to begin stim presentation at (motion wraps for the trial duration)
            % It will also increment on a PER EYE frame basis, but with a random starting point
            % Even though f0 and F could be the same thing, F starts at a random time point in the motion cycle, and indexes the motion
            % no matter how many frames are sampled (using mod(F,FramesFullCycle).
            % f0 only indexes the dot contrast values (for the temporal ramp) so must increment 1:PerEyeFrames
            F = round(rand*PerEyeFrames);
            
            while f0 <= PerEyeFrames % loop across stimulus frames (per eye)
                
                % Select left-eye image buffer for drawing:
                if useHardwareStereo
                    Screen('SelectStereoDrawBuffer', win, 0);
                end
                %%%%------------------------------------------------%%%%
                %               Draw left eye stimulus:
                %%%%------------------------------------------------%%%%
                
                % Draw dots:
                Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); % turn on alpha blending for lin superposition: see GarboriumDemo.m
                Screen('DrawTextures',win,DotsIdxL(:,mod(F,FramesFullCycle)+1), ...             % texture index
                    [], ...
                    dstRectsL(:,:,mod(F,FramesFullCycle)+1), ...    % select appropriate dot positions for current frame
                    [],[],TrialContrast(f0))                        % final argument is contrast
                
                % Superimpose the annulus:
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
                % Draw blue lines :
                Screen('FillRect', win, [0, 0, 1], blueRectLeftOn);
                Screen('FillRect', win, [0, 0, 0], blueRectLeftOff);
                
                % Select right-eye image buffer for drawing:
                if useHardwareStereo
                    Screen('SelectStereoDrawBuffer', win, 1);
                else
                    Screen('DrawingFinished', win);
                    [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); % update display on next refresh (& provide deadline)
                    
                    if missed > 0
                        missedFrames = missedFrames + 1;
                    end
                end
                
                %%%%------------------------------------------------%%%%
                %               Draw right eye stimulus:
                %%%%------------------------------------------------%%%%
                
                % Draw dots:
                Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); % turn on alpha blending for lin superposition: see GarboriumDemo.m
                Screen('DrawTextures',win,DotsIdxR(:,mod(F,FramesFullCycle)+1), ...             % texture index
                    [], ...
                    dstRectsR(:,:,mod(F,FramesFullCycle)+1), ...    % select appropriate dot positions for current frame
                    [],[], TrialContrast(f0))                        % final argument is contrast
                
                % Superimpose the annulus:
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
                % Draw blue lines:
                Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
                Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
                
                Screen('DrawingFinished', win);
                
                [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
                
                % keep record of any missed frames:
                if missed > 0
                    missedFrames = missedFrames + 1;
                end
                
                % Increment frame counters on a PER EYE basis:
                f0 = f0+1; % For contrast index
                F = F+1;   % For stimulus position index
                
            end % end of single trial stimulus loop
            
            % Keep fixation rings/lock on screen between the intervals for 1 s
            if Interval == 1
                Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
                Screen('DrawTexture',win, fixationRingTextureInner);
                Screen('DrawTexture',win, fixationRingTextureOuter);
                %Draw the fixation cross:
                Screen('FillRect',win,[0 0 0],fixationCross);
                Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
                Screen('BlendFunction', win, GL_ONE, GL_ZERO);
                Screen('DrawingFinished', win);
                Screen('Flip', win);
                WaitSecs(0.5);
            end
            
        end % End of loop across the 2 intervals.
        
        % Blank the screen while we await a response
        Screen('Flip', win);
        
        %%%%------------------------------------------------%%%%
        %               Check for response
        %%%%------------------------------------------------%%%%
        
        % Wait for a valid button press. If that response is correct, then indicate the response
        % as such. If not correct, response is coded as incorrect by default.
        % The loop will continue to wait until a valid button press is made
        % Responses can be made on the keyboard, or, if you're on the mac,
        % on the gamepad.
        
        % Wait for response if we haven't had one already ...
        validResp = false;
        % Continues to check for responses until  vald one is made.
        while ~validResp
            
            % Check for keyboard response
            [keyIsDown, ~, keyCode] = KbCheck();
            
            % Only use gamepad if you're on the mac:
            if ismac
                if ~keyIsDown % And only check if a keyboard response hasn't been made
                    [keyIsDown, ~, keyCode] = KbCheck(GamepadIndex); % Check for gamepad response
                end
            end
            
            if keyIsDown %If they've pressed a button on keyboard or gamepad,
                % determine if it is valid and proceed
                
                %If they respond 'interval 1' and the signal is in interval 1: **CORRECT**
                if any(keyCode(Interval1)) ...
                        && (trial_order(Trial, 1))      % Signal in interval 1
                    Responses(Trial) = 1;               % code response as correct
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(2000,0.4,0.05); Beeper(2000,0.4,0.05);
                    end
                    
                    % If they respond 'interval 2' and the signal is in interval 2: **CORRECT**
                elseif any(keyCode(Interval2)) ...
                        && ~(trial_order(Trial, 1))     % Signal not in interval 1
                    Responses(Trial) = 1;               % code response as correct
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(2000,0.4,0.05); Beeper(2000,0.4,0.05);
                    end
                    
                    % If they respond 'interval 1' and the signal is in interval 2: **INCORRECT**
                elseif any(keyCode(Interval1)) ...
                        && ~(trial_order(Trial, 1))     % Signal not in interval 1
                    % No need to recode response (coded as incorrect by default).
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(400,0.4,0.15); Beeper(400,0.4,0.5);
                    end
                    
                    % If they respond 'interval 2' and the signal is in interval 1: **INCORRECT**
                elseif any(keyCode(Interval2)) ...
                        && (trial_order(Trial, 1))      % Signal in interval 1
                    % No need to recode response (coded as incorrect by default).
                    validResp = true;
                    if GiveFeedback
                        % Play tones to indicate response made:
                        Beeper(400,0.4,0.15); Beeper(400,0.4,0.5);
                    end
                    
                elseif keyCode(RespQuit)
                    % break out of program if 'q' is pressed
                    ForcedQuit = true;
                    ExitGracefully(UsingVP, ForcedQuit)
                    
                end
            end
        end
        
        % write stimulus & response details to ResponseArray before updating staircase...
        % In the order the psi/trial, the stim level,
        % the 1st interval (0=noise, 1=signal) designation,
        % the response (1 or 0 for correct/wrong), the alphaHat value
        % and the error associated with that alphaHat value
        ResponseArray(Trial,:) = [trident, C, trial_order(Trial, 1), Responses(Trial), stair{trident}.alphaHat, stair{trident}.alphaHatSD];
        
        %---------------------------------------------------------%
        %                   Update the staircase:
        %---------------------------------------------------------%
        
        % update the staircase for the next trial, based on the subject's response
        
        % *** NOTE *** the staircases are updated at the end of the final trial AFTER the value
        % has been placed in 'ResponseArray', meaning that the very last value stored in the staircase struct arrays
        % (StairCase{trident}.alphaHat is the *TRUE* final threshold estimate, not the final value in the ResponseArray matrix
        
        
        stair{trident} = PsiUpdate( stair{trident}, Responses(Trial));
        stair{trident} = PsiStep( stair{trident} );
        
        
    end %end of loop across trials
    
    %%%%-----------------------------------------------------------%%%%
    %                       Save the results
    %%%%-----------------------------------------------------------%%%%
    
    % Store the data and some info about the run
    results.ResponseArray = ResponseArray;
    results.SubjectCode = SubjCode;
    results.RunNumber = runN;
    results.StairCase = stair;
    % For quick reference, set aside the frequency & amplitude parameters
    parameters.MotionFrequencyInHz = frequency;
    parameters.DisparityInArcMin = amplitude;
    
    % Save the data and parameters information!
    save(fName, 'parameters', 'results')
    
    missedFrames % print out number of missed frames during experiment
    
    % Finally, plot the staircases.
    figure
    Ylab = {'Coherence (%)'};
    Titles = {'Psi 1', 'Psi 2'};
    
    % All cont/coherence and alpha hat values *100 to convert to percentage
    % Remember, ResponseArray is in the order
    % the psi/trial, the stim level,
    % the 1st interval (0=noise, 1=signal),
    % the response (1 or 0 for correct/wrong), the alphaHat value,
    % and the alphaHat SDs.
    for X = 1:max(results.ResponseArray(:,1)) % equivalent to num_psi
        
        ix = (results.ResponseArray(:,1)==X);                                     % index for values of the current staircase
        ixCorr = (results.ResponseArray(:,1)==X & results.ResponseArray(:,4)==1); % index for correctly-answered trials for current staircase
        subplot(1,2,X)
        % Plot incorrect responses in white first:
        plot(results.ResponseArray(ix,2) .* 100, 'ko-', 'MarkerFaceColor', 'w', 'MarkerSize', 10)
        set(gca, 'box', 'off', 'TickDir', 'out', 'FontSize', 14)
        % Now determine correct responses and plot on top:
        corrResp = nan(length(results.ResponseArray),1);
        corrResp(ixCorr) = results.ResponseArray(ixCorr,2);
        hold on
        plot(corrResp(ix) .* 100, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
        % Now plot alpha hat values (threshold estimates) for each trial in current staircase:
        %plot(results.ResponseArray(ix,5) .* 100, 'r.-')
        % Also plot the Alpha Hat standard deviations:
        errorbar(results.ResponseArray(ix,5) .* 100, results.ResponseArray(ix,6) .* 100, 'r.-')
        xlabel('Trial')
        ylabel(Ylab)
        title(Titles{X})
        ylim([0 100])
        
    end
    suptitle(['Observer ', results.SubjectCode, ', ', ...
        parameters.Cue, ' ', parameters.DotColourType, ' Psi staircases, ', num2str(parameters.MotionFrequencyInHz), ...
        ' Hz, ', num2str(amplitude), ' arcmin, Run ', num2str(runN)])
    legend('Incorrect', 'Correct', 'alpha hat')
    
    % Close everything down & exit:
    ExitGracefully (UsingVP, ForcedQuit)
    
catch MException
    
    rethrow (MException)
    psychrethrow(psychlasterror)
    ExitGracefully (UsingVP, ForcedQuit)
    error('Error!')
    
end % End of try/catch statement
end % End of main function

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
function [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
    (Cue, trial_order, C, Trial, amplitude, frequency, AdjustXBy, PPD, Dots, texrect, imsize, parameters, TestScone)
disp(trial_order)
% This is where the pre-generated dot positions are loaded and the disparity added.
% (see GenerateIOVDdots_2.m & GenerateCDdots.m)
% These are picked at random (with replacement) from those already in the folder.
% Note that trial_order is entered as a scalar to define the order of the current interval, where 0 = noise, 1=signal

% Indicate how many existing sets of dot stimuli are in the subfolder 'stimuli' for loading later.
StimInFolder = 5;

% Define the max. change in disparity (horiz shift of dots) for the 2 cues:
% the numerator is in arcmin (but converted to pixels) akin to the amplitude of the sine wave
IOVD_disparity = amplitude/60 * PPD;
CD_disparity = amplitude/60 * PPD;

% Convert frequency into a string so we can load the appropriate stimulus:
IOVDstrfreq = num2str(frequency);
IOVDstrfreq(IOVDstrfreq == '.') = '_'; % replace decimal point with underscore
CDstrfreq = num2str(frequency);
CDstrfreq(CDstrfreq == '.') = '_';

% Draw a random value to load the stimuli with:
RandDraw = ceil(rand*StimInFolder);

%%%%-------------------------------------------------------------------------%%%%
%                     Shift dot trajectories: add disparity
%%%%-------------------------------------------------------------------------%%%%

% *** ---------------------------------- ***
%               Set up CD
% *** ---------------------------------- ***
if strcmp(Cue, 'CD')
    
    FileNameCD = fullfile('stimuli',['CD_dots_', CDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
    CD_stim = load(FileNameCD); % includes parameters file and dot positions matrices
    % Just cut the parameters back to a smaller struct:
    CDparams = CD_stim.DotParams;
    
    % *** set up dot coherence ***
    % Now we'll determine dot coherence according to the staircase
    NumDots = round(CDparams.NumDots * C); % modify number of coherent dots by value by current staircase setting
    
    % Now that we know the number of coherent dots, set up some random indices that will determine which dots are 'randomised'
    % this is necessary to ensure that dots of both colours are made noise dots equally.
    % We also want these indices to be consistent across frames in a trial, so they are set here.
    % Remember that we want to index the NOISE dots to change (ie=1), not the signal dots,
    % so it may seem a bit counterintuitive that noise =1, coherent =0
    NoiseDotIndices = Shuffle([true(CDparams.NumDots-NumDots,1); false(NumDots,1)]); %The 1s in here should sum to IOVDparams.NumDots-NumDots: the number of noise dots.
    
    % Define the number of frames in a full cycle:
    FramesFullCycle = CDparams.FramesFullCycle;
    
    % Assign the pre-generated dot position matrices.
    % First, we duplicate the left eye dot positions to make the right eye dot positions.
    % At this stage all dots are identical for the two eyes (disparity has not been added yet).
    dot_posL = CD_stim.dot_posL;
    dot_posR = CD_stim.dot_posL;
    
    % set aside the disparity
    CDparams.disparity = CD_disparity;
    % For reference also store the disparity in Arcmin
    CDparams.DisparityInArcmin = (CD_disparity/PPD)*60;
    
    % set up a few extra parameters of the sine waves for each cue:
    % the period, in sec
    CDparams.period = 1/CDparams.frequency;
    
    % the angular frequency of the sine waves:
    CDparams.angFreq = 2 * pi * CDparams.frequency;
    
    % Length of the sine wave:
    % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
    % end of the cycle and the beginning of the next cycle.
    % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
    CDparams.t = linspace(CDparams.period/CDparams.FramesFullCycle, CDparams.period, CDparams.FramesFullCycle);
    
    % Now make one full cycle of the sine wave, no matter the frequency
    % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
    % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
    CDparams.SineWave =  CDparams.disparity * sin(CDparams.angFreq * CDparams.t); %disparity * sin(angFreq * t)
    
    % assign dstRects: destination rect matrices
    dstRectsL = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
    dstRectsR = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
    
    % Shift dot trajectory: add disparity
    for fr = 1:CDparams.FramesFullCycle
        
        % determine dot coordinates: remember, y position does not change: horiz disparity only
        % Update dot position according to sinsoidal trajectory on each (per-eye) frame
        % Left eye: -sin
        dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ... % the size of the dot texture
            dot_posL(:,1,fr) - CDparams.SineWave(fr), ...  % the x positions
            dot_posL(:,2,fr))';                            % the y positions
        % Right eye: +sin
        dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
            dot_posR(:,1,fr) + CDparams.SineWave(fr), ...
            dot_posR(:,2,fr))';
        
        % Now, randomise the dot positions of the noise dots when coherence <1.0
        if (C < 1)
            %Left eye:
            dstRectsL(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(CDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ... % random x positions
                rand(CDparams.NumDots-NumDots,1) * imsize)';               % random y positions
            %right eye:
            dstRectsR(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(CDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ...
                rand(CDparams.NumDots-NumDots,1) * imsize)';
        end
        
    end
    
    % Now, if it's a NOISE interval, we simply suffle the frames temporally to give us the CD control stimulus.
    % This gives us a stimulus with the same distribution of disparities but no consistent/smooth
    % change in disparity over time
    % This should be all we need to change for the CD noise stimulus control
    if trial_order == 0
        ShuffledFrames = randperm(CDparams.FramesFullCycle);
        dstRectsL = dstRectsL(:,:,ShuffledFrames);
        dstRectsR = dstRectsR(:,:,ShuffledFrames);
    end
    % If it is a SIGNAL interval we don't need to change anything
    
    % Now we have the dot position indices (the dstRects), define the dot texture indices.
    % These are simply half white (Dots(1) and half black (Dots(2)) in each eye.
    
    %NOTE: since fixing the Right eye IOVD flicker bug, DotsIdx matrices are now 2D with num_dots * FramesFullCycle dimensions.
    % So we need to set up the matrix accordingly. We won't do it like the below any more:
    %DotsIdxL = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
    %DotsIdxR = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
    
    % Because in CD the dots only last one frame & then move, it doesn't really matter how their colour is assigned (as long as there is equivalent
    % numbers of black and white AND they are the same across the eyes.
    % So we can set them up in the same way as IOVD, except they are not later sorted by colour (using MinDist), and they are identical in the 2 eyes:
    
    % ALSO, here is the switch where we assign whether the dots are coloured (S-cone) or achromatic.
    % All spatial and temporal parameters for coloured/achromatic dots are the same.
    
    if ~TestScone % for the achromatic dots
        %Left eye:
        DotsIdxL = [repmat(Dots(1),1,CDparams.NumDots/2); repmat(Dots(2),1,CDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
        DotsIdxL = repmat(reshape(DotsIdxL,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
        %Right eye:
        DotsIdxR = [repmat(Dots(1),1,CDparams.NumDots/2); repmat(Dots(2),1,CDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
        DotsIdxR = repmat(reshape(DotsIdxR,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
        
    elseif TestScone % Or for the coloured dots.
        % Note that we measure isoluminance separately for the 2 eyes so each eye has its own pair of dot textures
        
        %Left eye:
        DotsIdxL = [repmat(Dots(3),1,CDparams.NumDots/2); repmat(Dots(4),1,CDparams.NumDots/2)]; %Dots(3) = purple; Dots(4)=yellow
        DotsIdxL = repmat(reshape(DotsIdxL,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
        %Right eye:
        DotsIdxR = [repmat(Dots(5),1,CDparams.NumDots/2); repmat(Dots(6),1,CDparams.NumDots/2)]; %Dots(5) = purple; Dots(6)=yellow
        DotsIdxR = repmat(reshape(DotsIdxR,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
        
    end
    
    % Set aside details of the stimulus parameters.
    % Most of these are already stored in the pre-generated stimulus files,
    % but it is probably worthwhile to have them here again in a single output file
    if trial_order == 1
        parameters.StimulusNumberEachTrial{Trial,1} = FileNameCD; % Store the filename of the randomly-drawn stimulus for SIGNAL
    elseif trial_order == 0
        parameters.StimulusNumberEachTrial{Trial,2} = FileNameCD; % Store the filename of the randomly-drawn stimulus for NOISE
    elseif trial_order == -1
        % -1 indicates that it is during set up, so set aside parameters
        parameters.CDparameters = CDparams;                       % Store the CD parameters once, during set-up
    end
    
    % *** ---------------------------------- ***
    %               Set up IOVD
    % *** ---------------------------------- ***
elseif strcmp(Cue, 'IOVD')
    
    FileNameIOVD = fullfile('stimuli',['IOVD_dots_colour_', IOVDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
    IOVD_stim = load(FileNameIOVD); % includes parameters file and dot positions matrices
    % Just cut the parameters back to a smaller struct:
    IOVDparams = IOVD_stim.DotParams;
    
    % *** set up dot coherence ***
    % Now we'll determine dot coherence according to the staircase
    NumDots = round(IOVDparams.NumDots * C); % modify number of coherent dots by value by current staircase setting
    
    % Now that we know the number of coherent dots, set up some random indices that will determine which dots are 'randomised'
    % this is necessary to ensure that dots of both colours are made noise dots equally.
    % We also want these indices to be consistent across frames in a trial, so they are set here.
    % Remember that we want to index the NOISE dots to change (ie=1), not the signal dots,
    % so it may seem a bit counterintuitive that noise =1, coherent =0
    NoiseDotIndices = Shuffle([true(IOVDparams.NumDots-NumDots,1); false(NumDots,1)]); %The 1s in here should sum to IOVDparams.NumDots-NumDots: the number of noise dots.
    
    % Define the number of frames in a full cycle:
    FramesFullCycle = IOVDparams.FramesFullCycle;
    
    % set aside the disparity
    IOVDparams.disparity = IOVD_disparity;
    % For reference also store the disparity in Arcmin
    IOVDparams.DisparityInArcmin = (IOVD_disparity/PPD)*60;
    
    % set up a few extra parameters of the sine waves for each cue:
    % the period, in sec
    IOVDparams.period = 1/IOVDparams.frequency;
    
    % the angular frequency of the sine waves:
    IOVDparams.angFreq = 2 * pi * IOVDparams.frequency;
    
    % Length of the sine wave:
    % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
    % end of the cycle and the beginning of the next cycle.
    % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
    IOVDparams.t = linspace(IOVDparams.period/IOVDparams.FramesFullCycle, IOVDparams.period, IOVDparams.FramesFullCycle);
    
    % Now make one full cycle of the sine wave, no matter the frequency
    % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
    % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
    IOVDparams.SineWave =  IOVDparams.disparity * sin(IOVDparams.angFreq * IOVDparams.t); %disparity * sin(angFreq * t)
    
    % assign dstRects: destination rect matrices
    dstRectsL = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
    dstRectsR = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
    
    % Shift dot trajectory: add disparity.
    % Remember that with IOVD the dot positions for each eye are NOT identical (before the horizontal shifts are added)
    % as they are with CD, that is why we have 2 matrices for the 2 eyes
    % Dot lifetime is determined by how many dot positions are re-randomised with each frame (see GenerateIOVDdots.m)
    
    % Determine a value to modulate the sine wave by.
    % In this way the sinusoidal disparity is either subtracted or added to the dot positions with each frame.
    % This way we can make the noise/control stimulus with a very simple switch
    % The two columns in SineModulate are for all dots in each eye: col 1 for LE, col 2 for RE.
    
    if trial_order == 0 % for a noise trial
        % Make vector of 1s & -1s
        SineModulate = ones(IOVDparams.NumDots,2);
        SineModulate(1:2:end,:) = -1; % this bit will make stimuli the IOVD control: both eyes have dots moving both directions
        
        %**** need random 1/2 to be black, 1/2 white; because color is assigned in alternate manner
        % motion direction cannot also.
        % If motion direction is assigned as above, in an alternate manner, the control has dots of opposite colours moving in
        % opposite directions. We can't have that.
        % So, re-assign SineModulate. Now it is a 2 * num_dots vector, upper half are -1s, lower half +1s
        % Remember that there is no relationship between dot position on the screen & its location in the matrix, so making the SineModulate
        % matrix in this way won't matter; the motion directions will still be assigned randomly in the IOVD control.
        % And because color is assigned in an alternate manner, that means equal numbers of black/white dots will
        % be going both monocular directions.
        
        SineModulate = ones(IOVDparams.NumDots,2);
        SineModulate(1:end/2,:) = -1;
        
        % *** Not doing the below any more ...
        %(Also, define indices to select a random sample of NumDots/2 from each eye.
        %All these dots will be presented to BOTH eyes to ensure the same amount of dots per eye.
        %These indices should stay constant across frames to ensure continuous motion of the dots (so they are assigned here).
        %Equal numbers of dots are sampled from each eye.
        %To do this, we simply re-assign the random sample of dots from the L eye to the R eye.
        %The resulting R eye matrix then has 50% of dots from the L & R eyes (because 50% of dots in the R eye matrix have not been re-assigned).
        %Then we simply make LeyeMatrix = ReyeMatrix so that the same pattern is displayed to both eyes.
        %LERandDotSample = Shuffle([true(IOVDparams.NumDots/2,1); false(IOVDparams.NumDots/2,1)]); )... ***
        
    elseif trial_order % when trial_order = 1 or -1 (or any value other than 0 really)
        % subtract left eye, add right eye for MID: signal interval
        SineModulate = ones(IOVDparams.NumDots,2);
        SineModulate(:,1) = -1; %-1 * sin for the left eye only in actual IOVD; +1 * sin for right eye
    end
    
    for fr = 1:IOVDparams.FramesFullCycle
        
        % determine dot coordinates: remember, y position does not change: horiz disparity only
        % Update dot position according to sinsoidal trajectory on each (per-eye) frame
        
        % Left eye: -sin
        dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ...                                      % the size of the dot texture
            IOVD_stim.DotPosBinoc{1}(:,1,fr) + SineModulate(:,1) * IOVDparams.SineWave(fr), ... % x positions
            IOVD_stim.DotPosBinoc{1}(:,2,fr))';                                                 % y positions
        
        %Right eye: +sin
        dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
            IOVD_stim.DotPosBinoc{2}(:,1,fr) + SineModulate(:,2) * IOVDparams.SineWave(fr), ...
            IOVD_stim.DotPosBinoc{2}(:,2,fr))';
        
        % *** Don't do the below any more:
        % AS usual we present 2 sets of dots to L/R eyes,
        % But the direction of motion is alternately left or right within a single eye, & this is determined by
        % 'SineModulate' parameter....
        
        % If it's a noise interval, we take a random sample of NumDots dots from both eyes,
        % these dots are displayed to BOTH eyes simultaneously
        %         if trial_order == 0
        %             %if a noise interval, L & R are identical
        %             dstRectsR(:,LERandDotSample,fr) = dstRectsL(:,LERandDotSample,fr); %reassign 50% of dots from L eye into R eye
        %             dstRectsL(:,:,fr) = dstRectsR(:,:,fr); %now make LE pattern = RE pattern
        %         end % ***
        
        %Now, randomise the dot positions of the noise dots when coherence <1.0
        %L & R eye matrices both exist by this point in time
        %The random indices for the NOISE dots were determined above, according to IOVDparams.NumDots-NumDots
        if (C < 1)
            
            %Left eye:
            dstRectsL(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(IOVDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ... % random x positions
                rand(IOVDparams.NumDots-NumDots,1) * imsize)';               % random y positions
            %right eye:
            dstRectsR(:,NoiseDotIndices,fr) = CenterRectOnPoint(texrect, ...
                rand(IOVDparams.NumDots-NumDots,1) * imsize + AdjustXBy, ...
                rand(IOVDparams.NumDots-NumDots,1) * imsize)';
        end
        
    end
    
    % Now we have the dot position indices (the dstRects), define the dot texture indices.
    % We use separate vectors for L & R eye just in case you want to split the eyes by colour.
    % *** This is where the dot color must be assigned in an alternate fashion for L & R eye
    % to make (mostly) sure that neighbouring dots across the eyes are of a different color
    % However you do it, colour must be assigned oppositely to the 2 eyes for this to work
    
    % Previously we only assigned dot color in alternating fashion for IOVD, not the control.
    % This was because we previously used the same set of dots in the 2 eyes for the control,
    % so if color was assigned in an alternating fashion then the dots would flicker between the 2 colours as the
    % dots were presented in alternating fashion to the 2 eyes.
    % Now, because we use a different set of dots, it is ok for the dots to be assigned in alternating manner
    % across the eyes. This also maintains the method used to reduce binocular correlations in the new IOVD algorithm, applying it to the control also.
    
    % To solve the R eye flicker bug associated with colours changing across dot lifetime, the DotIdx matrices
    % must now be a 2D matrix with FrmsFullCycle * num_dots dimensions
    % Hence, both DotsIdxL & DotsIdxR must now be indexed by the current frame in the cycle as the dots are drawn
    % (as was done with the dot position matrices).
    
    % Generate the alternating matrices of dot colours. Note that they need to alternate.
    % If we just assigned top half of the matrix to balck, bottom half white (for eg.), then that would mean
    % mostly white dots would go one direction in the IOVD control stimulus, mostly black in the other direction, for example.
    %(this is because dots are now positioned in alternating strips across the 2 eyes)
    % We can't have that; in the control we would want equal numbers of white/black dots going in either direction.
    
    % ****
    % ALSO, here is the switch where we assign whether the dots are coloured (S-cone) or achromatic.
    % All spatial and temporal parameters for coloured/achromatic dots are the same.
    % ****
    
    % For the achromatic dots:
    if ~TestScone
        
        %Left eye: WBWBWB ....
        DotsIdxL = [repmat(Dots(1),1,IOVDparams.NumDots/2); repmat(Dots(2),1,IOVDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
        DotsIdxL = repmat(reshape(DotsIdxL,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
        %Right eye: BWBWBW....
        DotsIdxR = [repmat(Dots(2),1,IOVDparams.NumDots/2); repmat(Dots(1),1,IOVDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
        DotsIdxR = repmat(reshape(DotsIdxR,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
        
        %Or to split by colour:
        %DotsIdxL = repmat(Dots(2),1,IOVDparams.NumDots); %white dots only
        %DotsIdxR = repmat(Dots(2),1,IOVDparams.NumDots); %black dots only
        
    elseif TestScone
        
        % Left eye: WBWBWB ....
        DotsIdxL = [repmat(Dots(3),1,IOVDparams.NumDots/2); repmat(Dots(4),1,IOVDparams.NumDots/2)]; % Dots(3) = purple; Dots(4)=yellow
        DotsIdxL = repmat(reshape(DotsIdxL,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
        % Right eye: BWBWBW....
        DotsIdxR = [repmat(Dots(5),1,IOVDparams.NumDots/2); repmat(Dots(6),1,IOVDparams.NumDots/2)]; %Dots(5) = purple; Dots(6)=yellow
        DotsIdxR = repmat(reshape(DotsIdxR,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
        
    end

    % Now sort these dot colour matrices according to the indices established earlier that account for dot lifetime
    % Because dot colour was assigned in an opposite manner above, we sort the matrices using the same indices.
    % This retains the opposite polarity assignment across eyes
    DotsIdxR = DotsIdxR(IOVD_stim.MinDist);
    DotsIdxL = DotsIdxL(IOVD_stim.MinDist);
    
    % Set aside details of the stimulus parameters.
    % Most of these are already stored in the pre-generated stimulus files,
    % but it is probably worthwhile to have them here again in a single output file
    if trial_order == 1
        parameters.StimulusNumberEachTrial{Trial,1} = FileNameIOVD; % Store the filename of the randomly-drawn stimulus for SIGNAL
    elseif trial_order == 0
        parameters.StimulusNumberEachTrial{Trial,2} = FileNameIOVD; % Store the filename of the randomly-drawn stimulus for NOISE
    elseif trial_order == -1
        % -1 indicates that it is during set up, so set aside parameters
        parameters.IOVDparameters = IOVDparams;                     % Store the CD parameters once, during set-up
    end
    
end

end