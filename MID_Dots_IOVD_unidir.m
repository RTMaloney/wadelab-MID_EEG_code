function MID_Dots_IOVD_unidir (dirn)

%This is a demo-type function that will loop through an IOVD stimulus repeatedly until
%a key is pressed.
%The inter-ocular velocity difference (IOVD) cue involves opposite lateral motion in the two eyes, using independent (uncorrelated)
%sets of moving dots. The dots within a single eye are correlated across time, so there is some coherent motion,
%but they are not correlated between the eyes, so there should be no changing disparity. This is an uncorrelated (or time-correlated) random dot stereogram.
%It will take some time in generating the stimuli (depending on factors such as
%dot density & the frequency of the sine wave modulation), so please be patient.
%This is because each dot needs to be separated by a minimum amount (so each dot must be checked against every other dot)
%Dots also have to be re-randomised every 3 frames to limit their dot lifetime, & further, this needs to be done twice, once for each eye.
%R Maloney July 2015
%
%   Version _2 : updated example, showing the new IOVD algorithm (see WorkOutNewIOVDdotSpacing.m and GenerateIOVDdots_2.m)
% R Maloney 4 Dec 2015
% Version _2 of this code was developed to remove as many spurious binocular matches from the 2 eyes' images as possible,
% giving us a 'purer' IOVD stimulus.
% Assign dots in L & R eye within different strips for IOVD.
% Dot positions are assigned in adjacent horizontal strips in the L & R eye to prevent overlap
% between the 2 eyes.
% The horizontal strips are 2 dot widths wide.
% Dot positions can be anywhere within this strip for each eye.
% Dot sigma (SD of the Gaussian envelope), and hence dot size/width, is now important for determining dot positions.
% After placing the dot positions within the strips, the smallest distances between dots is determined
% in order to deal with the dots that might fall close to the borders.
% The dots in the R eye are sorted according to their smallest distance to any dot in the L eye.
% Once a L eye dot has been determined to be closest to one R eye dot, it is removed from the
% distance computations that follow for the next R eye dot.
% That way a L eye dot is never assigned more than once, which could be possible if more than one
% L eye dot has their closest distance to a single R eye dot.
% After this sorting, the color (ie Black/white) of those dots must be assigned in alternating sequence:
% that way if a L eye dot is near the border of the strips and is close to a
% R eye dot position, then at least it will be assigned the opposite polarity
% (and hence be less likely to be paired binocularly). But the dot color assignment must happen
% prior to the presentation of the dots.

% Updated Jan 2018: displays a continuous
% UNIDIRECTIONAL sawtooth wave (towards or away according to dirn), that wraps at the end of
% the cycle

% Inputs: dirn = -1 for towards: negative sawtooth wave
%         dirn = +1 for away: positive sawtooth wave (the default).

% Set the flag for the direction of motion: away
if nargin < 1 || isempty (dirn)
    dirn = 1;
end

%Define some of the display parameters:
PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx (46.4 on the Propixx)

%If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = false;

%Switch for whether you want the annulus superimposed over the dots:
DrawAnnulus = true;

%And if you want to draw the fixation cross/rings:
DrawRings = true; 

%To independently switch off L or R eyes:
LeftEyeContrast = 1;
RightEyeContrast = 1;

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

if ~ismac
    jheapcl; %clear the java heap space.
end

% Define the dot texture, a square-shaped sheet of dots.
%Make the texture the same size as the height of the screen
imsize = screenRect(4) ;

%Define dot density in dots/(deg^2):
dot_dens_per_deg2 = 1; %2.273;

%compute number of dots in the total available area
num_dots = round(dot_dens_per_deg2 * (imsize/PPD)^2); %this many dots in the full dot field

%Just check whether num_dots is odd or even: (important for when contrast polarity is assigned)
%If odd, -1 to make it even
if mod(num_dots,2) %if odd, mod = 1
    num_dots = num_dots-1;
end

%specify dot size:
dot_sigma_in_degrees = 0.05; %0.1 for amblyopes; %size of SD of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; %sigma in pixels
dotsize = round(dot_sigma * 10); %make the dots some multiple of sigma
%NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
%It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.

%Define the minimum spacing between the PEAKS of the dot positions (in pixels):
SpatialJitter = round(0.5 * PPD); %+ dotsize/2;
%NOTE: add dot radius (ie dotsize/2) ensures the EDGES of the dots are separated by the minimum, but there may not be space in the matrix!

% *** Work out the dot strips. We will make the dot strips 2 dot widths wide.
% Positions can vary +/-1 dot width from the midpoint of each strip.
% Determine the dot width.
% At full width, half maximum (FWHM), there are approx 2.355 sigma under the Gaussian.
% So dot_sigma * 2.355 will give us the dot width at half the max height.
% Multiply by PPD to get result in pixels.
% Multiply by 2 to give us approx the full dot width

Dot_width = dot_sigma_in_degrees * 2.355 * PPD * 2; %The approx dot width.
StripWidth = 2*Dot_width; %So the strips are (approx) 2 dot widths wide.

%Compute the midpoints of the strips:
midpts = Dot_width:StripWidth:imsize;

%Take left and right eye strip midpoints:
%use a cell array because these vectors may have different lengths
y_pos_midpts{1} = midpts(1:2:end); %Left eye strip midpoints
y_pos_midpts{2} = midpts(2:2:end); %Right eye strip midpoints

%make the dot profile:
x = (1-dotsize)/2:(dotsize-1)/2;
%x = x/PPD; %rescale x into vis angle
[x,y] = meshgrid(x,x);
y = -y;
[a,r] = cart2pol(x,y);

%This gives us a white-peaked dot (+ve contrast polarity)
env = exp(-r.^2/(2*dot_sigma.^2));%gaussian window
env = env./max(max(abs(env))); % normalize peak to +/- 1
env2 = -env; %make a negative version to give us a black dot (-ve contrast polarity)

%set up the raised cosine annular window.
%specify parameters for the annulus:
inrad = PPD * 1;% inner radius of annulus (in pixels), for fixation spot
outrad = PPD * 12/2; %outer radius of annulus (in pixels)
% define extent of spatial raised cosine at edge of aperture (in pixels)
cos_smooth = dotsize; %make it one dot size wide
%This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2; %double the texture size
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

%%%%-------------------------------%%%%
%       Set up fixation
%%%%-------------------------------%%%%

%Set up the fixation cross or spot:
%This is drawn directly to the screen using Screen('FillRect')
crossWidth = 2;
crossHeight = 10;
fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

%Make the fixation lock ring:
% We have an inner one around fixation and an outer one right on the edge of screen.
% These could probably be defined as a single texture (rather than 2) but I think that will complicate matters with the alpha-blending settings.
% (they are complicated enough already)
ringRadiusInner = PPD*0.5;                % ring surrounding fixation
ringRadiusOuter = screenRect(4)/2;        % outer edge (radius) of the ring: the edge of the screen
ringWidthInner = ringRadiusInner - PPD/4; % 1/4 of a degree thick
ringWidthOuter = ringRadiusOuter - PPD/3; % 1/3 of a degree thick

%Make the ring. It's in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
%Define the ring:
xx = (1-imsize)/2:(imsize-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
% make the alpha mask for the rings, inner and outer.
ring_alphaInner = ((r>ringWidthInner+1) & (r<ringRadiusInner-1)); % Make the alpha mask a tiny bit thinner than the ring itself.
ring_alphaOuter = ((r>ringWidthOuter+1) & (r<ringRadiusOuter-1));

%Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

try %Start a try/catch statement, in case something goes awry with the PTB functions
    
    %----------------------------
    % Set up the screen
    %----------------------------
    
    % initialization of the display
    AssertOpenGL;
    % Open PTB onscreen window: We request a 32 bit per colour component
    % floating point framebuffer if it supports alpha-blending. Otherwise
    % the system shall fall back to a 16 bit per colour component framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    %Set the color range to be normalised between 0 and 1 (rather than 0-255):
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
    %Initialise the Vpixx device:
    
    if UsingVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        %The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        Datapixx('EnableVideoStereoBlueline');
        Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses
        
        if Datapixx('IsViewpixx3D') %If it's the Viewpixx3D
            
            %Datapixx('EnableVideoLcd3D60Hz');
            Datapixx('DisableVideoLcd3D60Hz'); %actually, disabling seems better for reducing crosstalk...

            subjectData.DisplayType = 'Viewpixx3D'; %set aside the device type for reference
            Datapixx('RegWr');
            
        elseif Datapixx('IsPropixx') %if it's the Propixx DLP projector
            
            subjectData.DisplayType = 'PROpixx'; %set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram',0); %set to normal RGB video processing for driving the LEDs & DLP MMDs
            %Datapixx('RegWr');
            
            %Modify the per-eye crosstalk on the PROpixx.
            %Apparently this cross-talk correction only works when using RB3D video mode,
            %where the red/blue channels contain the left/right eye greyscale images (which we are not using).
            %Datapixx('SetPropixx3DCrosstalkLR', 1);
            %Datapixx('SetPropixx3DCrosstalkRL', 1);
            Datapixx('RegWrRd'); %seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        end
    end
    %No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
    % Open an on screen (grey) window and configure the imaging pipeline
    %Info about the 'blueline' mechanism for synching to the 3D glasses:
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
    
    %Define the 'blue line' parameters
    blueRectLeftOn   = [0,                 windowRect(4)-1, windowRect(3)/4,   windowRect(4)];
    blueRectLeftOff  = [windowRect(3)/4,   windowRect(4)-1, windowRect(3),     windowRect(4)];
    blueRectRightOn  = [0,                 windowRect(4)-1, windowRect(3)*3/4, windowRect(4)];
    blueRectRightOff = [windowRect(3)*3/4, windowRect(4)-1, windowRect(3),     windowRect(4)];
    
    %HideCursor;
    %Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we shouldn't use
    %Screen(‘LoadNormalizedGamma’) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
    %The PROpixx device should have a linear lUT built in, but we will add this here for completeness.
    %R_gamma =
    %G_gamma =
    %B_gamma =
    %PsychColorCorrection('SetEncodingGamma', win, [1/R_gamma, 1/G_gamma, 1/B_gamma]);
    %raise priority level:
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
    %Query the screen refresh rate:
    ifi = Screen('GetFlipInterval',win); %in sec
    
    %Set the alpha-blending:
    %We want a linear superposition of the dots should they overlap:
    %Just like the Gabors in GarboriumDemo.m (see there for further info).
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);
    % We also want alpha-blending for smooth (anti-aliased) dots...
    %not sure how this will conflict with the above command
    %about the linear superposition of dots... but it doesn't seem to cause problems
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    %Make the textures for the dots, 1 for black and white
    Dots(1) = Screen('MakeTexture',win,env,[],[],2); %white dot
    Dots(2) = Screen('MakeTexture',win,env2,[],[],2); %black dot
    
    %Generate the annulus texture:
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
    
    %%%%-------------------------------------------------------------------------%%%%
    %       Set up values of the sinusoidal change in horizontal dot disparity
    %%%%-------------------------------------------------------------------------%%%%
    
    
    %Fixed parameters:
    %If the frequency of the sine wave is specified in Hz (ie cycles/sec)
    %then other units of time must also be in SECONDS!!!
    frequency = 4; %1.1 Hz based on RM & JA thresholds; temp frequency of motion, in Hz. We want 1/2 cycle in 500 ms, so 1 Hz
    period = 1/frequency; %in sec
    angFreq = 2 * pi * frequency; %for the sine wave (the periodic sine wave is always some multiple of 2*pi)
    TrialDuration = 0.5; %trial duration, IN SEC
    IOVD_disparity = 128/60 * PPD; %102 arcmin based on RM & JA thresholds; disparity in pixels, defined in arcmin (akin to the amplitude of the sine wave)
    %the effective frame rate PER EYE: since each eye is stimulated on successive video frames
    %Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
    PeyeFR = RefreshRate/2; %Per eye f.r. is always half the absolute f.r.
    %Frames MUST be a function of the screen RefreshRate, and the duration, otherwise the timing will be wrong.
    %frames = round(RefreshRate*TrialDuration); %also round to integer number of frames
    
    %Determined parameters:
    %The no. of frames for a full cycle:
    %Remember that each eye samples the motion vector at the same point, so this is in terms of per-eye frame rate
    %(otherwise the different eyes would sample successive points on the trajectory/sine wave: one eye would lead the other)
    FrmsFullCyc = round(PeyeFR*period); %must have integer no. of frames
    
    %Length of the sine wave:
    %The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
    %end of the cycle and the beginning of the next cycle.
    %So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
    t = linspace(period/FrmsFullCyc, period, FrmsFullCyc); %One complete cycle, IN SEC, in steps of per-eye frames
    
    %Now make one full cycle of the sine wave, no matter the frequency
    %Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
    %Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
    
    % *** This is where the single direction is determined *** RM, Jan 2018
    SineWave =  dirn * IOVD_disparity * sawtooth(angFreq * t); %IOVD_disparity * sin(angFreq * t)
    
    %assign dstRects matrix
    dstRectsL = zeros(4, num_dots,FrmsFullCyc);
    dstRectsR = zeros(4, num_dots,FrmsFullCyc);
    
    
    %%%%-------------------------------------------------------------------------%%%%
    %                   Set up initial dot positions
    %%%%-------------------------------------------------------------------------%%%%
    
    
    %determine random dot x and y positions, with minimum separation 'SpatialJitter'
    %ie pseudo-random positioning to prevent dot overlap/clustering.
    
    % Dot positions are assigned by first randomly selecting one of the strip midpoints for that eye.
    % Once the random midpoint is selected, the random value r is generated and added.
    % It is somewhere between -1 dot width & +1 dot widths, given by 'r = (2*rand-1) * Dot_width'
    % The dots are still checked & re-randomised if they fail the minimum separation 'spatial jitter'.
    
    disp('Setting up dot positions, please wait...')
    tic
    for EachEye = 1:2 %loop across both eyes, L then R
        
        CurrDot=1;
        dot_pos = zeros(num_dots,2); %assign dot location matrix
        y_pos = y_pos_midpts{EachEye}; %select the midpoints, left, then right
        
        while CurrDot <= num_dots
            
            if CurrDot == 1 %set random coordinates for very first dot
                r = (2*rand-1) * Dot_width; % *** gives random val betw -1 dot width & 1 dot width
                %Select a random strip midpoint for this eye, and randomly select 'y' position within the strip
                dot_pos(CurrDot,:) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %was imsize.*rand(1,2);
                CurrDot = CurrDot+1;
            else
                %set the next dot's random position
                r = (2*rand-1) * Dot_width; % *** gives random val betw -1 dot width & 1 dot width
                dot_pos(CurrDot,:) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %imsize.*rand(1,2);
                %find the smallest distance (in pixels) between the current dot and any other existing dot
                idx = 1:CurrDot-1; %index each existing dot except the current one
                d = min((dot_pos(idx,1)-dot_pos(CurrDot,1)).^2 + (dot_pos(idx,2)-dot_pos(CurrDot,2)).^2);
                d = sqrt(d);
                
                %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
                %This will continue until (at least) the minimum distance is met
                if d < SpatialJitter
                    r = (2*rand-1) * Dot_width;
                    dot_pos(CurrDot,:) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %imsize.*rand(1,2);
                else %if that minimum distance is met, move on to the next dot
                    CurrDot = CurrDot+1;
                end
            end
        end
        
        %Set aside dot positions for each eye:
        DotPosBinoc{EachEye} = dot_pos;
        
    end
    
    %%%%-------------------------------------------------------------------------%%%%
    %       replicate dot_pos by the number of frames needed (to allocate memory)
    %%%%-------------------------------------------------------------------------%%%%
    
    %Note that, previously, we were just swapping the x & y positions just generated for the L eye to produce
    %an essentially different set of positions for the R eye (& hence saving computation time), while still satisfying the 'spatial jitter' requirement.
    %We can't do that any more (because dots are now assigned in alternating strips for L & R).
    %The previous method did not guarantee the dots wouldn't overlap in L/R eyes, & they may have coincided along the diagonal anyway
    %DotPosBinoc{2} = repmat([dot_pos(:,2), dot_pos(:,1)],1,1,FrmsFullCyc); %for the right: swap x & y coordinates
    
    DotPosBinoc{1} = repmat(DotPosBinoc{1},1,1,FrmsFullCyc);
    DotPosBinoc{2} = repmat(DotPosBinoc{2},1,1,FrmsFullCyc);
    
    %%%%-------------------------------------------------------------------------%%%%
    %                   Determine dot lifetime
    %%%%-------------------------------------------------------------------------%%%%
    
    %We need to re-randomise 1/3 of all dots on every frame
    %This ensures all dots have a lifetime of 3 frames (or 50 ms, with a per-eye frame rate of 60 Hz).
    %So we set up the indices for (roughly) every 1/3 of the dots, leaving no dot unturned
    DotThirdIndices = round(linspace(1,num_dots,4));
    DotThirdIndices = [DotThirdIndices(1:3)', DotThirdIndices(2:4)'];
    DotThirdIndices(2:3,1) = DotThirdIndices(2:3,1)+1;
    
    %Repeat for the 2 eyes, since both eyes need different (uncorrelated) dots for IOVD
    for TwoEyes = 1:2
        
        dot_pos = DotPosBinoc{TwoEyes}; %re-assign dot_pos, for either the left (1) or right (2) eye
        y_pos = y_pos_midpts{TwoEyes}; %select the midpoints, left, then right
        
        %set dot positions for each frame in a full cycle
        for m = 1:FrmsFullCyc
            
            %Make the dots the same position as the previous frame.
            %1/3 will then be moved.
            if m~=1
                dot_pos(:,:,m) = dot_pos(:,:,m-1);
            end
            
            Curr3rd = mod(m,3)+1; %tells us whether this is a 1st, 2nd or 3rd frame, & determines which 3rd of dots to change
            CurrRows = DotThirdIndices(Curr3rd,:);
            dot_pos(CurrRows(1):CurrRows(2),:,m) = nan; %make the third we want to change NaNs to put them out of the calculations
            CurrDot = CurrRows(1);
            
            while CurrDot <= CurrRows(2) %go through all dots in this 3rd
                
                %set the next dot's random position
                r = (2*rand-1) * Dot_width; % *** gives random val betw -1 dot width & 1 dot width
                %Select a random strip midpoint for this eye, and randomly select 'y' position within the strip
                dot_pos(CurrDot,:,m) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r];                %find the smallest distance (in pixels) between the current dot and any other existing dot
                %Index all the existing dots except the current one to do this
                %This means excluding all the dots currently unassigned (defined as NaNs.
                CurrentEmpty = find(isnan(dot_pos(:,1,m)))'; %all rows currently empty
                idx = setdiff(1:num_dots, [CurrDot, CurrentEmpty]);
                
                %Find the smallest distance between the current dot and any other dot
                d = min((dot_pos(idx,1,m)-dot_pos(CurrDot,1,m)).^2 + (dot_pos(idx,2,m)-dot_pos(CurrDot,2,m)).^2);
                d = sqrt(d);
                
                %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
                %This will continue until (at least) the minimum distance is met
                if d < SpatialJitter
                    r = (2*rand-1) * Dot_width;
                    dot_pos(CurrDot,:,m) = [imsize.*rand, y_pos(ceil(rand*length(y_pos))) + r]; %imsize.*rand(1,2);
                else %if that minimum distance is met, move on to the next dot
                    CurrDot = CurrDot+1;
                end
                
            end %end of loop across the dots in the current 3rd
            
        end %end of loop across all frames in a full cycle
        
        %adjust the x positions so they are in the centre of the screen. Do this for all frames in one go:
        AdjustXBy = (screenRect(3) - screenRect(4))/2; %shift X dot positions by this amount to centre them, since image size <screenRect(3)
        dot_pos(:,1,:) = dot_pos(:,1,:) + AdjustXBy;
        
        DotPosBinoc{TwoEyes} = dot_pos; %set aside the newly-assigned dots for this eye
        
    end %end of loop across eyes
    
    %%%%-------------------------------------------------------------------------%%%%
    %                   Sort all dots by smallest distances between L/R eyes
    %%%%-------------------------------------------------------------------------%%%%
    
    %Remember this just sorts the order by which the dot positions are listed in the matrices so
    %that the L & R eye dots are paired according to their smallest distance to each other: this means we can
    %assign color to L & R eye in an alternate manner (eg BWBWBW...; WBWBW) meaning that any dots appearing close to one
    %another across the 2 eyes *should* have the opposite color/contrast polarity (*though there may still be the odd occasion
    %where this does not happen).
    
    % *** %
    % On any given frame 1/3 of dots change positions from their previous location.
    % This ensures there are the same amount of transient events across frames (even the 1st).
    % we also need to make sure these changes are tracked so a dot doesn't accidentally change colour during its lifetime
    % when dot colour is assigned.
    % Because of the sorting below, position in the MATRIX now matters when assigning dot colour, but 
    % the indices we retain in 'MinDist' allow us to keep track of that.
    % *** %

    %container for all min distance INDICES for each dot pair; for each frame
    % 1/3 of these indices change with each frame, ensuring each dot pairing lasts no longer than the 3 frame dot lifetime
    MinDist = nan(num_dots,FrmsFullCyc); 
    for F = 1:FrmsFullCyc %Need to do this for every frame
        
        Curr3rd = mod(F,3)+1; %tells us whether this is a 1st, 2nd or 3rd frame, & determines which 3rd of dots to change
        %NotCurr3rd = setdiff([1:3],Curr3rd)
        CurrRows = DotThirdIndices(Curr3rd,:);
        
        if F == 1
            %For the very first frame, we want to sort all of the dots
            dot_posL(:,:) = DotPosBinoc{1}(:,:,F);
            dot_posR(:,:) = DotPosBinoc{2}(:,:,F);
            D = pdist2(dot_posL(:,:), dot_posR(:,:));
            CurrRows = [1, num_dots];
        else
            % Make all values nan to take the 2/3 of dots we don't need to sort out
            dot_posL(:,:) = nan;
            dot_posR(:,:) = nan;
            %Then place the 1/3 of dot coordinates on this frame that we DO need to sort, back in
            dot_posL(CurrRows(1):CurrRows(2),:) = DotPosBinoc{1}(CurrRows(1):CurrRows(2),:,F);
            dot_posR(CurrRows(1):CurrRows(2),:) = DotPosBinoc{2}(CurrRows(1):CurrRows(2),:,F);
            %Now compute the distance matrix, excluding the 2/3 that do not need sorting
            D = pdist2(dot_posL(:,:), dot_posR(:,:));
            %Make MinDist the same as the previous frame (then we re-arrange 1/3)
            MinDist(:,F) = MinDist(:,F-1);
        end
        
        % Compute the distance matrix (unsorted):
        % Gives the distances between each pairwise combination of dots in each eye
        %D = pdist2(DotPosBinoc{1}(:,:,F), DotPosBinoc{2}(:,:,F));
        
        % Go through each ROW in the distance matrix, D, which correspond to each ROW in the L eye matrix, DotPosBinoc{1}
        % Once a dot in the R eye (DotPosBinoc{2}) has been assigned to a Left eye dot, it is removed from further distance computations.
        % So we would take the next smallest distance, then the next smallest, and so on.
        RE_ind_rep = true(num_dots,1); %set logical indices for RE dots in this third.
        
        for ix = CurrRows(1):CurrRows(2) %Loop across all dots in the current third
            
            %find the column with the smallest distance
            %MinDist indexes the COLUMNS of D and the ROWS of DotPosBinoc{2}
            %for each item 1..n of DotPosBinoc{1}
            MinDist(ix,F) = find(D(ix,:) == nanmin(D(ix,RE_ind_rep))); %Find the INDEX of the smallest distance between LE dot ix and all (unassigned/remaining) RE dots
            
            % Now update the RE index, to remove the dot with the smallest distance at this current (i'th)
            % iteration from further computations of the smallest distance for subsequent LE dots
            RE_ind_rep(MinDist(ix,F)) = false;
            
        end %end of loop across dots for current frame
        
        %So now we have the dots sorted according to their smallest possible distances.
        %(without replacement).
    
        % **** At this point, we previously would sort the dots in the dot position matrix according to the smallest distances.
        % However, we no longer need to do that. Because we record the indices of the dot pairings/sortings in MinDist,
        % we can just use that to sort the dot colour assignment matrices instead. Dot position matrices need not change.
        % Sorting the dot position matrix for the R eye only was part of the reason we encountered the dot colour flickering bug.        
        % *** So no longer doing this below *** :
        %DotPosBinoc{2}(CurrRows(1):CurrRows(2),:,F) = DotPosBinoc{2}(MinDist(CurrRows(1):CurrRows(2),F),:,F);
        
    end %end of sorting loop across frames
    toc
    %%%%-------------------------------------------------------------------------%%%%
    %                           Shift dot trajectory
    %%%%-------------------------------------------------------------------------%%%%
    
    %Now determine the lateral shift of the dots according to the sinusoidal modulation.
    %All dots for a single eye shift by the same amount on a single frame. This amount changes according to time (frame).
    %The shift in the left and right eyes is opposite in direction.
    
    %Determine and set aside dot positions for each frame for a whole cycle. This is done at the PER-EYE frame rate.
    %If we are presenting at less than 1 cycle some of these frames will never be used.
    %But if it is cycling through over multiple cycles (as in this demo), then they will all be used, & should be sampled from the sine wave & wrapped smoothly.
    %Lower frequencies will require more frames (& hence more time to compute them) because more frames(time) are needed to complete one full cycle.
    
    %Whether we present the IOVD control or actual IOVD can all be controlled by the variable
    %'SineModulate':
    
    % subtract left eye, add right eye for MID: signal interval
    SineModulate = ones(num_dots,2);
    SineModulate(:,1) = -1; %-1 * sin for the left eye only in actual IOVD; +1 * sin for right eye
    
    
    %SineWave = SineWave * 0; % do this to completely null the lateral motion
    
    for fr = 1:FrmsFullCyc
        
        %determine dot coordinates: remember, y position does not change: horiz disparity only
        %Update dot position according to sinsoidal trajectory or each (per-eye) frame
        %Left eye: -sin
        dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ... % the size of the dot texture
            DotPosBinoc{1}(:,1,fr) + SineModulate(:,1) * SineWave(fr), ...     % the x positions
            DotPosBinoc{1}(:,2,fr))';                      % the y positions
        %Right eye: +sin
        dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
            DotPosBinoc{2}(:,1,fr) + SineModulate(:,2) * SineWave(fr), ...
            DotPosBinoc{2}(:,2,fr))';
        
    end
    
    %Set up dot texture index for each location; half black, half white
    %We use separate vectors for L & R eye just in case you want to split the eyes by colour.
    % *** This is where the dot color must be assigned in an alternate fashion for L & R eye
    % to make (mostly) sure that neighbouring dots across the eyes are of a different color
    % However you do it, colour must be assigned oppositely to the 2 eyes for this to work 
    
    % To solve the R eye flicker bug associated with colours changing across dot lifetime, the DotIdx matrices
    % must now be a 2D matrix with FrmsFullCycle * num_dots dimensions
    % Hence, both DotsIdxL & DotsIdxR must now be indexed by the current frame in the cycle as the dots are drawn
    % (as was done with the dot position matrices).
    
    % Generate the alternating matrices of dot colours. Note that they need to alternate.
    % If we just assigned top half of the matrix to balck, bottom half white (for eg.), then that would mean 
    % mostly white dots would go one direction in the IOVD control stimulus, mostly black in the other direction, for example.
    %(this is because dots are now positioned in alternating strips across the 2 eyes)
    % We can't have that; in the control we would want equal numbers of white/black dots going in either direction.
    
    %Left eye: WBWBWB ....
    DotsIdxL = [repmat(Dots(1),1,num_dots/2); repmat(Dots(2),1,num_dots/2)]; %Dots(1) = white; Dots(2)=black
    DotsIdxL = repmat(reshape(DotsIdxL,1,num_dots)', 1, FrmsFullCyc); % Replicate for each frame
    %Right eye: BWBWBW....
    DotsIdxR = [repmat(Dots(2),1,num_dots/2); repmat(Dots(1),1,num_dots/2)]; %Dots(1) = white; Dots(2)=black
    DotsIdxR = repmat(reshape(DotsIdxR,1,num_dots)', 1, FrmsFullCyc); % Replicate for each frame
    
    %Or to split by colour:
    %DotsIdxL = repmat(Dots(1),1,num_dots); %white dots only
    %DotsIdxR = repmat(Dots(2),1,num_dots); %black dots only
    
    % Now sort these dot colour matrices according to the indices established earlier that account for dot lifetime
    % Because dot colour was assigned in an opposite manner above, we sort the matrices using the same indices.
    % This retains the opposite polarity assignment across eyes
    DotsIdxR = DotsIdxR(MinDist);  
    DotsIdxL = DotsIdxL(MinDist); 
    
    f = 0; %this value increases with each iteration, but on a per eye basis only
    %And it determines the phase of the motion at the given frame (ie place in the cycle)
    missedFrames = 0;
%     WaitSecs(2)
    KbCheck();
%     KbCheck();
%     WaitSecs(2)
    vbl = Screen('Flip',win); %sync vbl to start time
    while ~KbCheck
        
        % Select left-eye image buffer for drawing:
        if useHardwareStereo
            Screen('SelectStereoDrawBuffer', win, 0);
        end
        
        %%%%------------------------------------------------%%%%
        %               Draw left eye stimulus:
        %%%%------------------------------------------------%%%%
        
        %Draw dots:
        Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
        Screen('DrawTextures',win,DotsIdxL(:,mod(f,FrmsFullCyc)+1), ... % dot colour, now indexed each frame
            [],dstRectsL(:,:,mod(f,FrmsFullCyc)+1), ... % dot position
            [],[],LeftEyeContrast) %Final argument is contrast
        
        %Superimpose the annulus:
        if DrawAnnulus
            Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTexture',win,annulus);
            Screen('Blendfunction', win, GL_ONE, GL_ZERO);
        end
        
        %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
        if DrawRings
            
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTextureInner);
            Screen('DrawTexture',win, fixationRingTextureOuter);        %Draw the black fixation cross:
            Screen('FillRect',win,[0 0 0],fixationCross);
        end
        
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        
        %Draw blue lines:
        Screen('FillRect', win, [0, 0, 1], blueRectLeftOn);
        Screen('FillRect', win, [0, 0, 0], blueRectLeftOff);
        
        % Select right-eye image buffer for drawing:
        if useHardwareStereo
            Screen('SelectStereoDrawBuffer', win, 1);
        else %Not sure what would happen if this 'else' was actually true. I suspect something would go wrong, as stim would be presented twice.
            %But this is more or less how the Vpixx people do it in their 'DatapixxImagingStereoDemo'
            Screen('DrawingFinished', win);
            [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
            
            if missed > 0
                missedFrames = missedFrames + 1;
            end
        end
        
        %%%%------------------------------------------------%%%%
        %               Draw right eye stimulus:
        %%%%------------------------------------------------%%%%
        
        %Draw dots:
        Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
        Screen('DrawTextures',win,DotsIdxR(:,mod(f,FrmsFullCyc)+1), ... % dot colour, now indexed each frame
        [],dstRectsR(:,:,mod(f,FrmsFullCyc)+1), ... % dot position.
        [],[],RightEyeContrast) %Final argument is contrast
        
        %Superimpose the annulus:
        if DrawAnnulus
            Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTexture',win,annulus);
            Screen('Blendfunction', win, GL_ONE, GL_ZERO);
        end
        
        %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
        if DrawRings
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTextureInner);
            Screen('DrawTexture',win, fixationRingTextureOuter);        %Draw the fixation cross:
            Screen('FillRect',win,[0 0 0],fixationCross);
        end
        
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        
        %Draw blue lines:
        Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
        Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
        
        Screen('DrawingFinished', win);
        
        [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
        
        f = f+1 %increment counter for next frame
        
        %keep record of any missed frames:
        if missed > 0
            missedFrames = missedFrames + 1;
        end
        
    end
    f
    missedFrames
    
    %%%%------------------------------------------------%%%%
    %               Close down screen:
    %%%%------------------------------------------------%%%%
    %turn off the prioritisation:
    Priority( 0 ); %restore priority
    
    if UsingVP        % close down the ViewPixx or ProPixx
        Datapixx('DisableVideoScanningBacklight');
        if Datapixx('IsViewpixx3D')
            Datapixx('DisableVideoLcd3D60Hz');
        end
        Datapixx('RegWr');
        %Datapixx('Close'); %closing it here might cause it to crash?
    end
    
    %Close down the screen:
    Screen('CloseAll')
    
    Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing
    
    %Bring back the mouse cursor:
    ShowCursor();
    
catch MException
    
    if UsingVP        % close down the ViewPixx or ProPixx
        Datapixx('DisableVideoScanningBacklight');
        if Datapixx('IsViewpixx3D')
            Datapixx('DisableVideoLcd3D60Hz');
        end
        Datapixx('RegWr');
        %Datapixx('Close'); %closing it here might cause it to crash?
    end
    
    %Close down the screen:
    Screen('CloseAll')
    Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing
    
    rethrow (MException)
    psychrethrow(psychlasterror)
    
end %End of try/catch statement
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




