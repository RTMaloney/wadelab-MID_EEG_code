function MID_Dots_CD_unidir (dirn)

%This is a demo-type program that will loop continuously through a CD stimulus until a key is pressed.
%It shows how the changing disparity (CD) cue to motion in depth can be generated.
%This is done by drawing a dynamic random dot stereogram. In this, the dots across the 2 eyes are correlated (and only differ by the horiz. disparity)
%however they change position (and disparity) on each frame, so there is no continuous motion of any dot across time.
%Because there is a new set of dots required for each frame, drawing these dot positions & ensuring they are separated by the minimum distance 'SpatialJitter'
%may take some time. The left eye is generated first & then replicated for the right eye. At this stage they are identical. Then the disparity is added in different
%directions for the two eyes. Simply adding the disparity is less computationally intensive than drawing the initial dot positions.
%R Maloney July 2015

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
dot_sigma_in_degrees = 0.05; %0.075 for amblyopes %size of SD of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; %sigma in pixels
dotsize = round(dot_sigma * 10); %make the dots some multiple of sigma
%NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
%It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.

%Define the minimum spacing between the PEAKS of the dot positions (in pixels):
SpatialJitter = round(0.5 * PPD); %was 0.43* PPD  + dotsize/2;
%NOTE: add dot radius (ie dotsize/2) ensures the EDGES of the dots are separated by the minimum, but there may not be space in the matrix!

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

% Make the rings. Both are in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
% Define the ring:
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
            
            Datapixx('EnableVideoLcd3D60Hz');
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
    inrect = repmat(texrect', 1, num_dots);
    
    %%%%-------------------------------------------------------------------------%%%%
    %       Set up values of the sinusoidal change in horizontal dot disparity
    %%%%-------------------------------------------------------------------------%%%%
    
    %Fixed parameters:
    %If the frequency of the sine wave is specified in Hz (ie cycles/sec)
    %then other units of time must also be in SECONDS!!!
    frequency = 2; %1.4 Hz based on thresholds of RM and JA %temp frequency of motion, in Hz. We want 1/2 cycle in 500 ms, so 1 Hz
    period = 1/frequency; %in sec
    angFreq = 2 * pi * frequency; %for the sine wave (the periodic sine wave is always some multiple of 2*pi)
    TrialDuration = 0.5; %trial duration, IN SEC
    CD_disparity = 16/60 * PPD; %11 arcmin, based on RM and JA thresholds; disparity in pixels, defined in arcmin (akin to the amplitude of the sine wave)
    %the effective frame rate PER EYE: since each eye is stimulated on successive video frames
    %Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
    PeyeFR = RefreshRate/2; %Per eye f.r. is always half the absolute f.r.
    %Frames MUST be a function of the screen RefreshRate, and the duration, otherwise the timing will be wrong.
    frames = round(RefreshRate*TrialDuration); %also round to integer number of frames
    
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
    SineWave =  dirn * CD_disparity * sawtooth(angFreq * t); %CD_disparity * sin(angFreq * t)
    
    %assign dstRects matrices
    dstRectsL = zeros(4, num_dots,FrmsFullCyc);
    dstRectsR = zeros(4, num_dots,FrmsFullCyc);
    
    
    %%%%-------------------------------------------------------------------------%%%%
    %                   Sort out initial dot positions
    %%%%-------------------------------------------------------------------------%%%%
    
    
    %determine random dot x and y positions, with minimum separation 'SpatialJitter'
    %ie pseudo-random positioning to prevent dot overlap/clustering
    
    disp('Setting up dot positions, please wait...')
    
    %assign space to dot location matrices for all frames
    %just do the left eye for now; they are currently identical (disparity hasn't been added)
    dot_posL = zeros(num_dots,2,FrmsFullCyc);
    
    %set dot positions for each frame in a full cycle
    tic
    for m = 1:FrmsFullCyc
        
        CurrDot=1;
        dot_pos = zeros(num_dots,2); %assign dot location matrix
        
        while CurrDot <= num_dots
            
            if CurrDot == 1 %set random coordinates for very first dot
                dot_pos(CurrDot,:) = imsize.*rand(1,2);
                CurrDot = CurrDot+1;
            else
                %set the next dot's random position
                dot_pos(CurrDot,:) = imsize.*rand(1,2);
                %find the smallest distance (in pixels) between the current dot and any other existing dot
                idx = 1:CurrDot-1; %index each existing dot except the current one
                d = min((dot_pos(idx,1)-dot_pos(CurrDot,1)).^2 + (dot_pos(idx,2)-dot_pos(CurrDot,2)).^2);
                d = sqrt(d);
                
                %Now if that distance is smaller than the minimum permitted distance, re-randomise the dot coordinates
                %This will continue until (at least) the minimum distance is met
                if d < SpatialJitter
                    dot_pos(CurrDot,:) = imsize.*rand(1,2);
                else %if that minimum distance is met, move on to the next dot
                    CurrDot = CurrDot+1;
                end
            end
        end
        
        %Set aside the dot position matrix for each frame (left eye only)
        dot_posL(:,:,m) = dot_pos;
        
    end
    
    %adjust the x positions so they are in the centre of the screen. Do this for all frames in one go:
    AdjustXBy = (screenRect(3) - screenRect(4))/2; %shift X dot positions by this amount to centre them, since image size <screenRect(3)
    dot_posL(:,1,:) = dot_posL(:,1,:) + AdjustXBy;
    
    %Duplicate L eye to give the R eye:
    dot_posR = dot_posL;
    toc
    
    %%%%-------------------------------------------------------------------------%%%%
    %                     Shift dot trajectory: add disparity
    %%%%-------------------------------------------------------------------------%%%%
    
    %Now determine the lateral shift of the dots according to the sinusoidal modulation.
    %All dots for a single eye shift by the same amount on a single frame. This amount changes according to time (frame).
    %The shift in the left and right eyes is opposite in direction.
    
    %Determine and set aside dot positions for each frame for a whole cycle. This is done at the PER-EYE frame rate.
    %If we are presenting at less than 1 cycle some of these frames will never be used.
    %But if it is cycling through over multiple cycles (as in this demo), then they will all be used, & should be sampled from the sine wave & wrapped smoothly.
    %Lower frequencies will require more frames (& hence more time to compute them) because more frames(time) are needed to complete one full cycle.
    f = 0; %this value increases with each iteration, but on a per eye basis only
    for fr = 1:FrmsFullCyc
        
        %determine dot coordinates: remember, y position does not change: horiz disparity only
        %Update dot position according to sinsoidal trajectory on each (per-eye) frame
        dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, dot_posL(:,1,fr) - SineWave(mod(f,FrmsFullCyc)+1), dot_posL(:,2,fr))';
        dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, dot_posR(:,1,fr) + SineWave(mod(f,FrmsFullCyc)+1), dot_posR(:,2,fr))';
        f = f+1;
        
    end
    
    %Set up dot texture index for each location; half black, half white
    %We use separate vectors for L & R eye just in case you want to split the eyes by colour
    DotsIdxL = [repmat(Dots(1),1,num_dots/2), repmat(Dots(2),1,num_dots/2)];
    DotsIdxR = [repmat(Dots(1),1,num_dots/2), repmat(Dots(2),1,num_dots/2)];
    %Or to split by colour:
    %DotsIdxL = repmat(Dots(1),1,num_dots); %white dots only
    %DotsIdxR = repmat(Dots(2),1,num_dots); %black dots only
    
    f = 0; %this value increases with each iteration, but on a per eye basis only
    %And it determines the phase of the motion at the given frame (ie place in the cycle)
    missedFrames = 0;
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
        Screen('DrawTextures',win,DotsIdxL,[],dstRectsL(:,:,mod(f,FrmsFullCyc)+1),[],[],1)
        
        %Superimpose the annulus:
        if DrawAnnulus
            Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTexture',win,annulus);
            Screen('Blendfunction', win, GL_ONE, GL_ZERO);
        end
        
        %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
        Screen('DrawTexture',win, fixationRingTextureInner);
        Screen('DrawTexture',win, fixationRingTextureOuter);
        %Draw the black fixation cross:
        Screen('FillRect',win,[0 0 0],fixationCross);
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
        Screen('DrawTextures',win,DotsIdxR,[],dstRectsR(:,:,mod(f,FrmsFullCyc)+1),[],[],1)
        
        %Superimpose the annulus:
        if DrawAnnulus
            Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
            Screen('DrawTexture',win,annulus);
            Screen('Blendfunction', win, GL_ONE, GL_ZERO);
        end
        
        %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
        Screen('DrawTexture',win, fixationRingTextureInner);
        Screen('DrawTexture',win, fixationRingTextureOuter);
        %Draw the fixation cross:
        Screen('FillRect',win,[0 0 0],fixationCross);
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        %Draw blue lines:
        Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
        Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
        
        Screen('DrawingFinished', win);
        
        [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
        
        f = f+1; %increment counter for next frame
        
        %keep record of any missed frames:
        if missed > 0
            missedFrames = missedFrames + 1;
        end
        
    end
    
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


