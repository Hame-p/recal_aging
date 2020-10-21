%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% spatial ventriloquist &  Recalibration paradigm
% 5 positions per modality
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;close all;clc;
%==========================================================================
% General parameters
%==========================================================================
% screen, sound, 
ARG.Do_trigger = 0; % set to 0 when testing on a different computer !
ARG.SNDCARD =  'Creative SBZ Series ASIO';
ARG.screenNr = 2;
ARG.Do_Eye = 0;
ARG.log_dir = 'C:/Users/eeglab/Desktop/CNSLAB/HP02/log';
ARG.paradigm = 'VRecal';


% timing
ARG.Timing_ITI = [900,300]; % ms Inter-trial timing [fixed delay, random component
ARG.Timing_Fix = [700,400]; % ms Interval between fixation onset and stimulus
ARG.Delay_Stim = [400,300]; % ms Delay between stimuli offset and response bar onset.

% stimulus specific
ARG.dist = 143; % distance between subject and screen in cm
ARG.VisRel= 2; % SD OF VISUAL DOTS
ARG.SND_scale = 0.2;
ARG.VIS_scale = 8; % distance between targets in degrees 
ARG.StimDur = 0.05;

ARG.tpcond = 11; % NUMBER OF REPEATS PER CONDITION 
ARG.totblk = 4;

Cond_name = {'AV','A','V'};
textstr = {'S', 'S', 'V'};

%==========================================================================
% Subject details
%==========================================================================
Subj = input('Subject: ','s');
ARG.Subj = Subj;
B = input(sprintf('Block (1-%d):', ARG.totblk));
ARG.Block = B;
EXIT_KEY = [81];
c = clock;
sname = sprintf('%s/VRecal_%s_B%d_%02d%02d-%02d%02d.mat',ARG.log_dir,Subj,B,c(3),c(2),c(4),c(5));
% EYETRACKER FILE NAME - needs to be short
etname = sprintf('VR%s%d.edf',Subj,B)
GetSecs;
%==========================================================================
% CREATE SCREEN
%==========================================================================
ssze(1) = 1; % screen width in cm
ssze(2) = 91; % screen height in cm
IMG_BGND = 60; % background luminance for gray screen
INIT_screen;
ARG.ifi = ifi;
fix_cord = [center-3 center+3]; % fixation dot

% screen size [SX, SY]
%=======================================================================
% initialize trigger port for EEG
%=======================================================================
INIT_trigger;

%=======================================================================
% initialize Eye tracker
%=======================================================================
dummymode=0;       % set to 1 to initialize in dummymode
INIT_Eye;

%=======================================================================
% visual stimuli
%=======================================================================

ARG.dots.ndots       = 200; % number of dots
ARG.dots.dot_w       = 0.12;  % width of dot (deg)
%ARG.dots.ppd = round((1/2)*SY*(2*ARG.dist*tan((2*pi/360)))/ssze(2));%54;    % pixels per degree % NOT CORRECT YET
ARG.dots.ppd = 28; % set to match with speakers

ARG.VisPos = (-2:2)*ARG.VIS_scale; % first is far-right, last is far-left.
ARG.locs = 1:5;

ppos = [nan ARG.VisPos];

%=======================================================================
% initialize sound presentation - here we need the 5 channel audio
%=======================================================================
ARG.Rate = 44100;
INIT_snd5;

% tone 
tone = sin([1:round(ARG.Rate*ARG.StimDur)]/ARG.Rate * 2*pi *1300);
% white noise
%tone = randn(1,round(ARG.Rate*ARG.StimDur));


% sounds are nchan x time (5 x time points)
% positions coded as 1 to 5 from right to left
[tone]=cosine_ramp(tone, 0.005, ARG.Rate);
for s=1:5    % must beware of the order since recalMEG is coded as 1:left, 5:right 
    SND{s}= zeros(5,length(tone));
    SND{s}(s,:) = tone;
end
% cf. Our speakers: SND{1} is right, SND{3} is center SND{5} is left
sswap = 5:-1:1; % swap sound positions from left-right.

%==========================================================================
% prepare response bar
%==========================================================================
% Make a bar of the perfect! size.
base_bar = [0 0 SX round(0.02*SY)];

% Put the bars in perfect! positions.
pres_bar = CenterRectOnPointd(base_bar, center(1), center(2));

% Set the bars to perfect! colors. First is for the bar, second for visual and
% third for auditory, last one is for the highlight pointer.
cols = [65 65 65; 
      0  0  255; 
      0  0  255; 
      255  0  0; 
      255  255  0]./255; % gray / blue / blue / red / yellow

% set cursor rectangle
curDimPix = round(0.03*SY/2); % length of line.
curLinePix = 6; % weight of the cursor line.
curCoords = [0 0 0 0 ; 0 0 -curDimPix curDimPix];

%==========================================================================
% conditions
%==========================================================================
% unisensory, multisenory, all congruencies
% modality: 0) bimodal 1) auditory, 2) visual
cond_file = sprintf('%s_S%s_B%d_cond.mat', ARG.paradigm, num2str(ARG.Subj, '%02d'), ARG.Block);
[condMAT, basecond] = create_condmat(ARG, ARG.log_dir, cond_file);
save(sname,'ARG','condMAT');
% condMAT: Ntrial x [vispos, soundpos, vispos-soundpos, modality]
% 1) vispos-soundpos is NaN for visual only trials
% 2) vispos-soundpos for sound only trials is the same as the previous AV
% trial
% 3) order is already randomized.

% number of trials per block:
Ntrial = size(condMAT, 1);
%==========================================================================
% Instructions
%==========================================================================

fprintf('Finished loading data\n');

DrawFormattedText(w, 'Antworten Sie mit der Maus', 'center', 300, [250 250 250]);
DrawFormattedText(w, 'Klick fuer Start', 'center', 500, [250 250 250]);
vbl = Screen('Flip', w);
if ARG.Do_Eye, Eyelink('StartRecording');  pause(0.05); end

HideCursor(w);
ShowCursor('CrossHair',w)
SetMouse(center(1), center(2),w);
Screen('ConstrainCursor', w, 1,[0 center(2) SX center(2)])
[clicks,x,y,whichButton] = GetClicks(w);


%==========================================================================
% main loop - trials
%==========================================================================
%% Set response matrix
response = nan(Ntrial, 12);

% set matrix for collecting stimuli onset and offset times
% columns: 1 - fixation onset, 2 - stimuli onset (visual, auditory), 3 -
% response bar onset
start_times = nan(Ntrial, 3);

fprintf('BLOCK: %d/%d \n', ARG.Block, ARG.totblk)
tic
for trial = 1:Ntrial
    
    Curr_stim = condMAT(trial,:); % vpos, apos, v-a, modality (0,1,2)
    trial_code = trial +1;
    if trial_code > 250
      trial_code = trial_code - 249;
    end
    HideCursor(ARG.screenNr)
    
    fprintf('trial %d of %d.  Mod: %s  Vp~Ap: %d  PosV: %d  PosA: %d  \n', ... 
      trial,Ntrial,Cond_name{Curr_stim(4)+1},ARG.VIS_scale*Curr_stim(3), ppos(Curr_stim(1)+1),ppos(Curr_stim(2)+1));
    
    % prepare sound 1) A, 2) V, 3)AV
    % for V only we use zero intensity
    if Curr_stim(4)==2
      PsychPortAudio('FillBuffer', pahandle, 0*SND{1}); % channelsxtime
    else
      PsychPortAudio('FillBuffer', pahandle, SND{sswap(Curr_stim(2))}); % channelsxtime
    end
    PsychPortAudio('Start', pahandle, 1, inf, 0);
    
    
    % prepare visual stimulus------------------------------------------
    if Curr_stim(4)~=1
      xydots = randn(ARG.dots.ndots,2)*ARG.VisRel*ARG.dots.ppd;
      xydots(:,1) = xydots(:,1)+ARG.VisPos(Curr_stim(1))*ARG.dots.ppd;
    else
      xydots = zeros(ARG.dots.ndots,2);
    end
    
    % wait until next trial start
    delay = (ARG.Timing_ITI(1)+rand*ARG.Timing_ITI(2))/1000;
    pause(delay);
    
    % present fixation dot
    Screen('FillOval', w, uint8(white), fix_cord);
    vbl=Screen('Flip', w);
    start_times(trial, 1) = vbl;
    
    % send trigger and store timing ----------------------------
    if ARG.Do_Eye, Eyelink('message',sprintf('TRIALID %d',trial)); end
    if ARG.Do_trigger, io64(cogent.io.ioObj,address,trial_code); end;
    t1 = GetSecs;
    pause(0.08);
    if ARG.Do_trigger, io64(cogent.io.ioObj,address,0); end;
    SetMouse(center(1), center(2),w);
    
    
    % wait period until targets
    delay = (ARG.Timing_Fix(1)+rand*ARG.Timing_Fix(2))/1000; % s delay
    pause(delay);
    
    % Stimuli ------------------------------------------------------
    
    StartTime = PsychPortAudio('RescheduleStart', pahandle, 0, 1);
    if ARG.Do_Eye, Eyelink('message', 'DISPLAY ON');  end
    if ARG.Do_trigger, io64(cogent.io.ioObj,address,1); end;
    
    Screen('FillOval', w, uint8(white), fix_cord);
    vbl=Screen('Flip', w);
    Screen('FillOval', w, uint8(white), fix_cord);
    Screen('DrawDots', w, transpose(xydots), ARG.dots.dot_w * ARG.dots.ppd, white, center,1);
    vbl=Screen('Flip', w);
    Screen('FillOval', w, uint8(white), fix_cord);
    Screen('Flip', w,vbl+ifi*3);
    start_times(trial, 2) = vbl;
    % wait period until response bar
    delay = (ARG.Delay_Stim(1)+rand*ARG.Delay_Stim(2))/1000; % s delay
    pause(delay);

    
    % ------------------------------------------------------------
    % wait for response with mouse
    % draw line
    % wait for click
    
    % set mouse cursor to random position on the bar.
    % SetMouse(round(rand * SX), center(2), w);
    % set mouse cursor to the middle of the bar.
    SetMouse(center(1), center(2), w);
    
    
    t2 = GetSecs; % this is when the subjects can start responding.
    start_times(trial, 3) = t2;
    % Present response bars and a mouse position.
    % continue to show the moving cursor until subject has responded to both stimuli.
    
    % set mouse cursor to random position on the bar.
    %         SetMouse(round(rand * SX), center, w);
    % set mouse cursor to the middle of the bar.
    SetMouse(center(1), center(2), w);
    
    butts = 0;
    while ~any(butts)
      [x, y, butts] = GetMouse(w);
      x = min(x, SX);
      x = max(x, 0);
      y = min(y, SY);
      % in the mean time show bar at mouse cursor.
      % Draw the bar: gray
      Screen('FillRect', w, cols(1, :), pres_bar);
      
      % Text output of mouse position draw in the centre of the screen
      DrawFormattedText(w, textstr{Curr_stim(4) + 1}, x - round(0.005*SX), ...
        center(2) - round(0.05*SY), cols(Curr_stim(4) + 2, :));
      
      % Draw a cursor line where the mouse cursor is to guide the
      % subject.
      Screen('DrawLines', w, curCoords,...
        curLinePix, cols(5, :), [x center(2)], 2);
      Screen('Flip', w);
    end
    
    %   [clicks,x,y,whichButton] = GetClicks(w);
    t3 =GetSecs;
    vbl=Screen('Flip', w);
    RT = t3-t2;
    
    [o1,o1,o2,stoptime] = PsychPortAudio('Stop', pahandle, 0);
    if ARG.Do_Eye, Eyelink('message', 'DISPLAY OFF');  end
    if ARG.Do_trigger, io64(cogent.io.ioObj,address,0); end;
        
    response(trial,:) = [Curr_stim,x,y,RT,whichButton,t1,t2,t3,0];
    save(sname,'ARG','response','start_times','basecond','-append');
end
elt = toc;
fprintf('it took %d minutes for this block.', elt/60)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DrawFormattedText(w, 'Ende', 'center', 500, [250 250 250]);
ShowCursor(w)
Screen('Flip', w);
pause(1);
Priority(0);
Screen('CloseAll');
PsychPortAudio('close')
if ARG.Do_Eye,
  Eyelink('message', 'ENDBUTTON');
  Eyelink('stoprecording');
  Eyelink('message', 'TRIAL OK');
  Eyelink('CloseFile');
  Eyelink('shutdown');
end
