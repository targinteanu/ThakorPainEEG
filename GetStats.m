%% Set Filepaths
clear 

home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)
postproDir = uigetdir;

cd(postproDir);
scanfiles = dir; 
[~,scanfiles] = listdlg_selectWrapper({scanfiles.name},'multiple');

cd(home); addpath(postproDir);
svloc = [postproDir,'/Stats ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];

%% select what to plot
meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);
plotOpts = {'Peak Freq (Hz)', 'Band Mean (\muV^2 s^2)', 'Band Density'};
fcnOpts = {@(w,P,band) peakFreq(w,P,band), meanAmp, ampDensity};
plotSel = listdlg_selectWrapper(plotOpts, 'single', 'Plot What?');
if sum(plotSel == 1:3)
    [bnd,bndname] = pickFrequency();
    yname = [bndname,' ',plotOpts{plotSel}];
    fcn0 = fcnOpts{plotSel};
    if sum(plotSel == [2,3])
        bnd = band2freqs(bnd, BandTableHz);
        fcn = @(Spec, EEG) frqFcnEpoch(Spec, EEG, @(w,P) fcn0(w,P,bnd));
        if plotSel == 3
            % Density 
            ylims = [0 1];
        else
            ylims = [];
        end
    elseif plotSel == 1
        % peak freq
        fcn = @(Spec, EEG) peakFreqEpoch(Spec, EEG, bnd, BandTableHz);
        if isa(bnd,'char') | isa(bnd,'string')
            ylims = band2freqs(bnd, BandTableHz);
        else
            ylims = bnd;
        end
    end
end
[~,testVars] = listdlg_selectWrapper(...
    {'TempStim', 'PinPrick', 'Pressure', 'BaselineOpen', 'BaselineClosed', 'BaselineIce'}, ...
    'multiple', 'Specify Variable');

%% construct plots for each subject 
table2baseline = @(tbl) [tbl.BaselineOpen('before experiment'),...
                         tbl.BaselineOpen('after experiment')];
p = .001; % uncertainty level 
timeBetweenEvents = 5; timeAfterLast = 30; % seconds 
baselineTransparency = .1;

dataTables = cell(length(scanfiles),2);
BLs = cell(length(scanfiles),1);
for subj = 1:length(scanfiles)
    fn = scanfiles(subj).name
    load(fn);

    % get stats on baseline 
    disp('calculating baseline statistics')
    BL_Epoch = table2baseline(Epoch_table);
    BL_EpochSpec = table2baseline(EpochSpec_table);
    BL = cell(size(BL_Epoch));
    for idx = 1:length(BL(:))
        if (~isempty(BL_Epoch{idx})) & (~isempty(BL_EpochSpec{idx}))
            [t,Y] = fcn(BL_EpochSpec{idx}, BL_Epoch{idx});
        end
        BL{idx} = cat(t,Y,3);
    end
    BL = cell2mat(reshape(BL,[],1)); BLY = BL(:,:,2);
    dof = size(BL,1) - 1; % deg of freedom 
    SE = std(BLY) / sqrt(size(BL,1)); % standard error
    BL_CI = mean(BLY) + [-1; 0; 1] .* tinv(p, dof)*SE;
    BLs{subj} = BL_CI;
    clear BL_Epoch BL_EpochSpec t Y idx dof SE BLY BL BL_CI

    % run calculations on desired variables 
    disp('calculating desired variables')
    [tY_table, EEG_table] = ...
        fcnTbl(EEG_table, Epoch_table, EpochSpec_table, fcn, testVars);

    dataTables{subj,1} = tY_table; dataTables{subj,2} = EEG_table;
    clear tY_table EEG_table Epoch_table EpochSpec_table Spec_table
end

% determine max trial duration 
disp('determining trial durations')
maxTrialDur = 0; % s
for subj = 1:size(dataTables,1)
    EEG_table = dataTables{subj,2}; 
    for r = 1:height(EEG_table)
        for c = 1:width(EEG_table)
            EEG = EEG_table{r,c}{1};
            init_time = [EEG.event.init_time];
            timeDiff = diff(init_time);
            sameTrial = [0, timeDiff <= timeBetweenEvents, 0]; 
                % ^ event i = same trial as event i-1
            trialBnd = diff(sameTrial); % 1=start, -1=end
            dur = init_time(trialBnd == -1) - init_time(trialBnd == 1);
            if isempty(dur)
                dur = 0;
            end
            maxTrialDur = max(maxTrialDur, max(dur));
            clear EEG init_time timeDiff sameTrial trialBnd dur
        end
    end
    clear EEG_table
end
maxTrialDur = maxTrialDur + timeAfterLast;

% segment into trials 
trialTables = cell(size(dataTables));
for subj = 1:size(dataTables,1)
    fn = scanfiles(subj).name
    EEG_table = dataTables{subj,2}; EEG_trial = EEG_table;
    tY_table  = dataTables{subj,1}; tY_trial  = tY_table;
    disp('segmenting into trials')
    for r = 1:height(EEG_table)
        for c = 1:width(EEG_table)
            EEG = EEG_table{r,c}{1};
            tY = tY_table{r,c}{1}; Y = tY(:,:,2); t = tY(:,:,1);
            init_time = [EEG.event.init_time];
            timeDiff = diff(init_time);
            sameTrial = [0, timeDiff <= timeBetweenEvents]; 
                % ^ event i = same trial as event i-1
            startEvs = EEG.event(~sameTrial);
            EEGs = repmat(EEG, size(startEvs));
            tYs = cell(size(startEvs));
            for idx = 1:length(startEvs)
                startT = startEvs(idx).latency/EEG.srate;

                disp(['Segment ',num2str(idx),' of ',num2str(length(startEvs)),...
                    ' (',num2str(100*idx/length(startEvs),3),'%)'])
                curEEG = pop_select(EEG, 'time', startT+[0,maxTrialDur]);
                    curEEG.xmin = curEEG.xmin + startT;
                    curEEG.xmax = curEEG.xmax + startT;
                    curEEG.times = curEEG.times + startT*1000;
                EEGs(idx) = curEEG;

                ti = (t >= startT) & (t < startT+maxTrialDur);
                curY = Y; curY(~ti) = nan; 
                curT = t; curT(~ti) = nan;
                tYs{idx} = cat(curT,curY,3);

                clear startT curEEG ti curY curT
            end
            EEG_trial{r,c} = {EEGs};
            tY_trial{r,c} = {tYs};
            clear EEG tY init_time timeDiff sameTrial EEGs tYs
        end
    end
    trialTables{subj,1} = tY_trial; trialTables{subj,2} = EEG_trial;
    clear EEG_table tY_table EEG_trial tY_trial
end

% organize for plotting 
for v = testVars
    disp(['setting up plots for ',v{:}]);

    disp('counting all trials')
    trialCount = zeros(size(trialTables,1),1);
    for subj = 1:size(trialTables,1)
        EEG_trial = trialTables{subj,2}; EEG_trial = makeSubtbl(EEG_trial, v);
        for r = 1:height(EEG_trial)
            trialCount(subj) = length(EEG_trial{r,1}{:});
        end
    end
    H = sum(trialCount > 0); trialTables = trialTables((trialCount > 0), :);
    W = max(trialCount);

    disp('plotting figures')
    figure; sgtitle([yname,' ',v{:}]);
    for subj = 1:H
        fn = scanfiles(subj).name
        strend = find(fn == '-'); strend = find(diff(strend)==1); strend = max(1, strend(1));
        ylbl = fn(1:strend);

        BL = BLs{subj};
        EEG_trial = trialTables{subj,2}; EEG_trial = makeSubtbl(EEG_trial, v);
        tY_trial  = trialTables{subj,1}; tY_trial  = makeSubtbl(tY_trial,  v);
        idx1 = 1;
        for r = 1:height(EEG_trial)
            EEGs = EEG_trial{r,1}{:};
            tYs = tY_trial{r,1}{:};
            for trl = 1:length(EEGs)
                EEG = EEGs(trl); chlocs = EEG.chanlocs;
                tY = tYs{trl}; Y = tY(:,:,2); t = tY(:,:,1);
                chanSig = sum( (Y < BL(:,1)) | (Y > BL(:,3)) );

                idx2 = (subj-1)*W + idx1;
                subplot(H,W,idx2);
                ttl = [EEG_trial.Properties.RowNames{r},' trial ',num2str(trl)];
                plotEvents(EEG);
                title(ttl); ylim(ylims); ylabel(ylbl); xlabel('time (s)');

                for chan = 1:size(Y,2)
                    if chanSig(chan)
                        tc = t(:,chan); Yc = Y(:,chan); 
                        plot(tc, Yc, ...
                            'Color',chanColor(chlocs(chan),chlocs),...
                            'LineWidth',1);
                    end
                end
                legend({chlocs.labels})
                for chan = 1:size(Y,2)
                    if chanSig(chan)
                        tc = t(:,chan); Yc = Y(:,chan); 
                        fill([min(tc),min(tc),max(tc),max(tc)],...
                             [BL(chan,1),BL(chan,3),BL(chan,3),BL(chan,1)],...
                             chanColor(chlocs(chan), chlocs),...
                             'FaceAlpha',baselineTransparency);
                    end
                end

                idx1 = idx1 + 1;
            end
        end
    end
    clear idx1 indx2 EEG_trial tY_trial EEGs tYs EEG tY t Y ttl ylbl strend ...
          chanSig chan chlocs tc Yc;
end

%% helper functions 

function subtbl = makeSubtbl(tbl, vars)
    subtbl = tbl(:, ismember(tbl.Properties.VariableNames, vars));
end

function rgb = chanColor(chloc, chlocs)
    allXYZ = [chlocs.X; chlocs.Y; chlocs.Z]';
    chXYZ  = [chloc.X;  chloc.Y;  chloc.Z ]';
    chXYZ = chXYZ - min(allXYZ);
    allXYZ = allXYZ - min(allXYZ);
    rgb = chXYZ./max(allXYZ);
end

function [bandrange, bandname] = pickFrequency()
    global BandTableHz
    freqOpts = [BandTableHz.Properties.RowNames; 'custom'];
    bandrange = listdlg_selectWrapper(freqOpts, 'single', 'Specify Frequency Band');
    if bandrange == length(freqOpts)
        % custom 
        bandrange = inputdlg({'Minimum Frequency (Hz):', 'Maximum Frequency (Hz):'},...
            'Specify Custom Frequency Band:');
        bandname = [bandrange{1},'-',bandrange{2},'Hz'];
        bandrange = arrayfun(@(i) str2double(bandrange{i}), 1:length(bandrange));
    else
        bandrange = freqOpts{bandrange}; bandname = ['\',bandrange];
    end
end

function [sel, listOut] = listdlg_selectWrapper(list, SelectionMode, PromptString)
    if nargin < 3
        PromptString = [];
        if nargin < 2
            SelectionMode = 'multiple';
        end
    end

    [sel, ok] = listdlg('ListString',list, 'SelectionMode',SelectionMode, 'PromptString',PromptString);
    while ~ok
        if strcmp(SelectionMode,'multiple')
            sel = questdlg('select all?');
            ok = strcmp(sel, 'Yes');
            if ~ok
                [sel, ok] = listdlg('ListString',list, 'SelectionMode',SelectionMode, 'PromptString',PromptString);
            else
                sel = 1:length(list);
            end
        else
            [sel, ok] = listdlg('ListString',list, 'SelectionMode',SelectionMode, 'PromptString',PromptString);
        end
    end
    listOut = list(sel);
end

%% key functions 

function [tblOut, toTestTbl, eegTbl, fig] = ...
    testTbl(toTestTbl, baseline, eegTbl, sttl, vars, rows)
    if nargin > 4
        toTestTbl = makeSubtbl(toTestTbl, vars, rows);
        eegTbl = makeSubtbl(eegTbl, vars, rows);
    elseif nargin < 4
        sttl = '';
        %{
        if nargin < 5
            yname = '';
            if nargin < 4
                ybound = [];
            end
        end
        %}
    end
    baseline = baseline(:,:,2);

    tblOut = toTestTbl; 
    fig = figure; sgtitle(sttl); 
    W = height(toTestTbl); H = width(toTestTbl); idx = 1;
    for c = 1:H
        for r = 1:W
            if ~isempty(toTestTbl{r,c}{1})
                curVar = toTestTbl{r,c}{1}(:,:,2);
                curTime = toTestTbl{r,c}{1}(:,:,1);
                curEEG = eegTbl{r,c}{1}(1);
                curTestOut = zeros(size(curVar));
                for chan = 1:size(curVar,2)
                    for t = 1:size(curVar,1)
                        [~,p] = ttest(baseline(:,chan),curVar(t,chan));
                        curTestOut(t,chan) = p;
                    end
                end
                tblOut{r,c} = {cat(3,curTime,curTestOut)};

                subplot(H,W,idx);
                ttl = [toTestTbl.Properties.VariableNames{c},' ',toTestTbl.Properties.RowNames{r}];
                plotWithEvents(curTime, curTestOut, curEEG, [0,1], ttl, 'p');

            end
            idx = idx + 1;
        end
    end
end

function [tblOut, eegTbl, epochTbl, epochSpecTbl] = ...
    fcnTbl(eegTbl, epochTbl, epochSpecTbl, fcn, vars)
    if nargin > 4
        eegTbl       = makeSubtbl(eegTbl,       vars);
        epochTbl     = makeSubtbl(epochTbl,     vars);
        epochSpecTbl = makeSubtbl(epochSpecTbl, vars);
    end
    
    tblOut = epochTbl; 
    W = height(epochTbl); H = width(epochTbl); 
    for c = 1:H
        for r = 1:W
            curSpec = epochSpecTbl{r,c}{:};
            curEpoc = epochTbl{r,c}{:};
            if ~isempty(curEpoc)
                [t,Y] = fcn(curSpec, curEpoc);
                tblOut{r,c} = {cat(3,t,Y)};
                eegTbl{r,c} = {eegTbl{r,c}(1)};
            end
        end
    end
end


function plotEvents(EEG)
    event = EEG.event; srate = EEG.srate;
    hold on; 
    lbl_sw = false;
    for ev = event
        if ~isempty(ev.latency) & ~strcmp(ev.type,'boundary')
            initTime = ev.latency/srate;
            plot(initTime+[0,0], ybound, 'r', 'LineWidth',1.5);
            text(initTime, ybound(lbl_sw+1), ev.type);
            lbl_sw = ~lbl_sw;
        end
    end
end

function [times,Y] = frqFcnEpoch(epoch_Spec, epoch_EEG, fcn)
    Y = cell2mat( arrayfun(@(eegSpec) ...
        fcn(eegSpec.frequency1side, eegSpec.powerSpectrum), ...
        epoch_Spec(2:end), 'UniformOutput', false) );
    times = arrayfun(@(eeg) mean([eeg.xmin, eeg.xmax]), epoch_EEG(2:end));
    times = repmat(times,size(Y,1),1);
    times = times'; Y = Y';
end

function [times,Freq] = peakFreqEpoch(epoch_Spec, epoch_EEG, bnd, tbl)
% inspect peakFreq across epochs (time)
    % index 1 = original 
    nchan = epoch_Spec(1).nbchan;
    PAFs = zeros(nchan, length(epoch_Spec)-1, 2);
    for chan = 1:nchan
        PAFs(chan,:,1) = arrayfun(@(eegSpec) ...
            peakFreq(eegSpec.frequency1side, eegSpec.powerSpectrum(chan,:), ...
            bnd, tbl), ...
            epoch_Spec(2:end));
        PAFs(chan,:,2) = arrayfun(@(eeg) ...
            mean([eeg.xmin, eeg.xmax]), epoch_EEG(2:end));
    end
    times = PAFs(:,:,2)'; Freq = PAFs(:,:,1)';
end