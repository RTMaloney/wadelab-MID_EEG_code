function Results = Process_MID_SSVEP_bh (SubjCode, DD_MON_YYYY, varargin)

% Inputs: SubjCode, subject code as a string, eg. 'RM'
% DD_MON_YYYY: the date string as in the file name, eg '01_Jan_2018'
% varargin: Enter run numbers consecutively as separate variables, depending on which ones you
% want to explore: eg. enter 1,2,3 for runs 1-3, enter 2 for run 2 only, etc.
% If more than 1 run is parsed, it will be concatenated and an overall score provided.
% Remember, the conditions are coded in this way:
% 1. CD towards
% 2. CD away
% 3. IOVD towards
% 4. IOVD away

% We will plot the subplot figures the same size each time, so obtain the screen dimensions to determine the size/location.
% We can make them the full width of the screen and half its height.
ScrSz = get(0,'ScreenSize');

%%
% load the data and store in cell array 'dat'.
% One for each run parsed.
for ii = 1:length(varargin)
    % build the file name
    fName = fullfile('data', [SubjCode, '_MID_SSVEP_', DD_MON_YYYY, '_', num2str(varargin{ii}), '.mat']);
    dat{ii} = load (fName);
end

% Button press data is encoded in per-eye FRAMES.
% Create a vector of timepoints at which the button presses are sampled, in SEC
RespSamplePts = linspace(0, dat{1}.parameters.RunDurationSec, length(dat{1}.subjectData.Responses));

% Loop across each run, concatenate the button presses if more than 1 run.
% Set up containers for each of the relevant variables.
StimOnsets = [];
ButtonPressTimes = [];
ButtonPresses = [];
ConditionOrder = [];
% We concatenate each subsequent run (after the first) by adding 300 sec (or the run duration) to the timestamps.
% That is like having one big long continuous run, rather than several separate short runs that return to 0.
for ii = 1:length(varargin)
    
    % Extract stim onset times for each run:
    StimOnsets = [StimOnsets, ...
        RespSamplePts(~isnan(dat{ii}.subjectData.Responses(:,1))) ...
        + (ii-1) * dat{ii}.parameters.RunDurationSec]; % add run duration in sec to concatenate
    % simply index the timestamps wherever a stimulus onset is indicated
    % NOTE: we already have a vector of stimulus onset asynchronies in parameters.cSOAs,
    % but that also includes the onsets of the BLANK periods, which we don't need here.
    
    % Do the same for the times where a button is pressed.
    % Note that a button press will be sampled across a few frames. We will take the very first time stamp of this
    % response as the response onset.
    % Only accept as response if it occurs within 1 sec of stim OFFSET: so that's 1sec + stim duration
    ButtonPressTimes = [ButtonPressTimes, RespSamplePts(~isnan(dat{ii}.subjectData.Responses(:,2))) ...
        + (ii-1) * dat{ii}.parameters.RunDurationSec]; % add run duration in sec to concatenate
    
    % Set aside the actual responses. -1 = towards, +1 = away.
    ButtonPresses = [ButtonPresses; dat{ii}.subjectData.Responses((~isnan(dat{ii}.subjectData.Responses(:,2))),2)];
    
    % Pull out the order of conditions.
    % We take all the values that are not '5s'. 5 indicates a blank.
    ConditionOrder = [ConditionOrder, dat{ii}.parameters.ConditionOrder( ...
        dat{ii}.parameters.ConditionOrder<5)];
    
end

%%
% Now we simply go through each of the 4 types of stim probe and check the response (if any).
% Define a valid window to accept responses within. We will say 1 sec + the probe duration (probably quite generous).
ValidResponseWindow = 1 + dat{1}.parameters.ProbeLengthSec;
for x=1:4
    
    indx = find(ConditionOrder==x); % Pull out indices for each condition type.
    NumStim(x) = length(indx); %Set aside the total number for that condition: should be the same for each.
    
    % Look across each probe.
    rts = nan(1,NumStim(x)); % Container for reaction times
    responses = zeros(1,NumStim(x)); % For responses. Label each as wrong, and change if correct.
    % We will keep a tally of missed responses: because a 'wrong' direction response
    % is not quite the same as a non-response (missed response), when we treat them as a 2AFC task.
    NumMissedResp(x) = 0;
    for y=1:NumStim(x)
        
        Btn_idx = find(ButtonPressTimes >= StimOnsets(indx(y)) ...         % If a button is pressed following the stim onset
            & ButtonPressTimes < StimOnsets(indx(y))+ValidResponseWindow); %...but before the end of the response window
        
        if ~isempty(Btn_idx) % If a response was made for this stimulus...
            rts(y) = ButtonPressTimes(Btn_idx(1)) - StimOnsets(indx(y)); %set aside the reaction time
            % And indicate when the response was correct:
            if ButtonPresses(Btn_idx(1)) == -1  % a 'towards' response
                if x == 1 || x == 3
                    responses(y) = 1; % mark as correct response
                end
            elseif ButtonPresses(Btn_idx(1)) == 1 % an 'away' response
                if x == 2 || x == 4
                    responses(y) = 1; % mark as correct response
                end
            end  
        else % If there was no response made, increment the number of missed 
             % responses for this condition. 
            NumMissedResp(x) = NumMissedResp(x) + 1; 
            % And change the response to a nan, to take it out of the later calculations:
            responses(y) = nan;
        end
    end
    % Set aside the data for each condition.
    Results.ReactionTimes(:,x) = rts;
    Results.PropCorr(:,x) = responses;
    Results.NumStim = NumStim;
    Results.NumMissedResponses = NumMissedResp; 
    Results.NumResponses(x) = NumStim(x) - NumMissedResp(x); % Determine the number of actual RESPONSES to stimuli made (so n stim - n missed responses)
end

%% --- Calculate sensitivity ---
% Let's compute sensitivity for the 2 cues, in the form of d'
% The equation adds the z-scores of the towards and away proportions correct
% for CD and IOVD, and divides by sqrt(2).

% Note that we cannot compute d prime if P(C) = 1.0. (or zero, but hopefully that will never happen!)
% So in that event, we use P(C) = (N trials-1)/N
Accu = nanmean(Results.PropCorr)
% Note that here we use the number of valid responses, rather than the total number of stimuli.
% Hopefully there will be very little difference in these numbers...
Accu(Accu==1) = (Results.NumResponses(Accu==1) - 1) / Results.NumResponses(Accu==1);

d_prime(1) = (norminv(Accu(1)) + norminv(Accu(2))) / sqrt(2); % For CD 
d_prime(2) = (norminv(Accu(3)) + norminv(Accu(4))) / sqrt(2); % For IOVD

% Compute asymptotic variance estimates for d' (see Sewell & Smith, 2011; Gourevitch and Galanter, 1967).
% Note that again, here we use the number of valid responses, rather than the total number of stimuli.
var_est(1) = Accu(1) * (1-Accu(1)) ...
    / ...
    (2 * Results.NumResponses(1) * normcdf(norminv(Accu(1)))^2) ...
    + ...
    Accu(2) * (1-Accu(2)) ...
    / ...
    (2 * Results.NumResponses(2) * normcdf(norminv(Accu(2)))^2); % For CD

var_est(2) = Accu(3) * (1-Accu(3)) ...
    / ...
    (2 * Results.NumResponses(3) * normcdf(norminv(Accu(3)))^2) ...
    + ...
    Accu(4) * (1-Accu(4)) ...
    / ...
    (2 * Results.NumResponses(4) * normcdf(norminv(Accu(4)))^2); % For IOVD

% Compute binomial standard errors for the proportions correct.
bin_SE =  sqrt((nanmean(Results.PropCorr) .* (1-nanmean(Results.PropCorr))) ./ Results.NumResponses);

% Set aside the sensitivities and the variance estimates
Results.Sensitivities = d_prime;
Results.Sensitivity_var_est = var_est;

%%
% Now, plot the data.
h = figure;
setpixelposition(h, [10, ScrSz(4)/3, ScrSz(3), ScrSz(4)/2])
suptitle([SubjCode, ', ', DD_MON_YYYY(DD_MON_YYYY~='_'), ', Run: ', num2str(cat(2,varargin{:}))])
% Proportions correct +/- binomial SES first:
subplot(1,3,1)
bar(Accu, 0.8, 'w')
hold on
errorbar(Accu, bin_SE, 'k', 'LineStyle', 'none')
ylim([0 1.05])
set(gca,'XTickLabel',{'CD Towards','CD Away','IOVD Towards', 'IOVD Away'}, ...
    'FontSize', 14, 'TickDir', 'out', 'box', 'off', 'XTickLabelRotation', 25)
ylabel('Proportion correct responses')

% Plot sensitivities (d')
subplot(1,3,2) 
bar(d_prime, 0.8, 'w')
hold on
errorbar(d_prime, var_est, 'k', 'LineStyle', 'none')
set(gca,'XTickLabel',{'CD','IOVD'}, ...
    'FontSize', 14, 'TickDir', 'out', 'box', 'off')
ylabel('Sensitivity (d'')')

% Plot RTs with standard error of the mean as error bars.
subplot(1,3,3) 
bar(nanmean(Results.ReactionTimes)*1000, 0.8, 'w')
hold on
errorbar(nanmean(Results.ReactionTimes)*1000, (nanstd(Results.ReactionTimes)*1000)./sqrt(NumStim), 'k', 'LineStyle', 'none')
set(gca,'XTickLabel',{'CD Towards','CD Away','IOVD Towards', 'IOVD Away'}, ...
    'FontSize', 14, 'TickDir', 'out', 'box', 'off', 'XTickLabelRotation', 25)
ylabel('Approx. reaction time (ms)')

