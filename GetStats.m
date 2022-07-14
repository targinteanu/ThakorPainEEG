%% Start eeglab
clear
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Set Filepaths
home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)
load("BrainwaveFrequencyTable.mat");
global BandTableHz
postproDir = uigetdir;

cd(postproDir);
scanfiles = dir; 
scanfiles = scanfiles((~strcmp({scanfiles.name},'.')) & (~strcmp({scanfiles.name},'..')));
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

%% baseline and calculations  
table2baseline = @(tbl) [tbl.BaselineOpen('before experiment'),...
                         tbl.BaselineOpen('after experiment')];
p = .05; % uncertainty level 
timeBetweenEvents = 5; timeAfterLast = 30; % seconds 
baselineTransparency = .1;

dataTables = cell(length(scanfiles),2);
BLs = cell(length(scanfiles),2);
for subj = 1:length(scanfiles)
    fn = scanfiles{subj}
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
        BL{idx} = cat(3,t,Y);
    end
    BL = cell2mat(reshape(BL,[],1)); BLY = BL(:,:,2);
    dof = size(BL,1) - 1; % deg of freedom 
    SE = std(BLY) / sqrt(size(BL,1)); % standard error
    BL_CI = mean(BLY) + [-1; 0; 1] .* tinv(p, dof)*SE;
    BLs{subj,1} = BL_CI; BLs{subj,2} = BL;
    clear BL_Epoch BL_EpochSpec t Y idx dof SE BLY BL BL_CI

    % run calculations on desired variables 
    disp('calculating desired variables')
    [tY_table, EEG_table] = ...
        fcnTbl(EEG_table, Epoch_table, EpochSpec_table, fcn, testVars);

    dataTables{subj,1} = tY_table; dataTables{subj,2} = EEG_table;
    clear tY_table EEG_table Epoch_table EpochSpec_table Spec_table
end

%% graph baselines of each subj 
fig(2) = figure; fig(1) = figure;
for subj = 1:size(BLs,1)
    fn = scanfiles{subj};
    strend = find(fn == '-'); strend = strend((diff(strend)==1)); strend = max(1, strend(1)-2);
    pname = fn(1:strend);
    pname = pname((end-2):end);

    % replace with more robust
    %EEG = dataTables{subj,2}{1,1}{1}(1); 
    EEG_table = dataTables{subj,2}; EEG = [];
    for r = 1:height(EEG_table)
        for c = 1:width(EEG_table)
            EEGrc = EEG_table{r,c}{1};
            if isempty(EEG)
                if ~isempty(EEGrc)
                    EEG = EEGrc;
                end
            end
        end
    end

    H = floor(sqrt(size(BLs,1)));
    W = ceil(size(BLs,1)/H);

    figure(fig(1));
    subplot(W,H,subj); 
    BL = BLs{subj,1};
    bar(BL(2,:));
    hold on; errorbar(BL(2,:), BL(3,:) - BL(2,:), ...
        '.k', 'LineWidth',1);
    xticks(1:length(BL(2,:)));
    xticklabels({EEG.chanlocs.labels}); xlabel('Channel'); ylabel(yname);
    title(['Subject ',pname,' Baseline']);
    ylim(ylims);

    figure(fig(2)); sgtitle(['Baseline ',yname]);
    subplot(W,H,subj); 
    topoplot(BL(2,:), EEG.chanlocs, 'maplimits', ylims, 'electrodes','labels'); colorbar;
    title(['Subject ',pname]);
end
clear fig BL EEG W H

%% determine max trial duration 
disp('determining trial durations')
maxTrialDur = 0; % s
for subj = 1:size(dataTables,1)
    EEG_table = dataTables{subj,2}; 
    for r = 1:height(EEG_table)
        for c = 1:width(EEG_table)
            EEG = EEG_table{r,c}{1};
            if ~isempty(EEG)
                evs = EEG.event; 
                evs = evs(arrayfun(@(ev) ~sum(strcmp(ev.type, {'-1','15','14','12','13'})), evs));
                init_time = [evs.init_time];
                timeDiff = diff(init_time);
                sameTrial = [0, timeDiff <= timeBetweenEvents, 0];
                % ^ event i = same trial as event i-1
                trialBnd = diff(sameTrial); % 1=start, -1=end
                dur = init_time(trialBnd == -1) - init_time(trialBnd == 1);
                if isempty(dur)
                    dur = 0;
                end
                maxTrialDur = max(maxTrialDur, max(dur));
            end
            clear EEG init_time timeDiff sameTrial trialBnd dur evs
        end
    end
    clear EEG_table
end
maxTrialDur = maxTrialDur + timeAfterLast;

%% segment into trials 
trialTables = cell(size(dataTables));
for subj = 1:size(dataTables,1)
    fn = scanfiles{subj}
    EEG_table = dataTables{subj,2}; EEG_trial = EEG_table;
    tY_table  = dataTables{subj,1}; tY_trial  = tY_table;
    disp('segmenting into trials')
    for r = 1:height(EEG_table)
        for c = 1:width(EEG_table)
            EEG = EEG_table{r,c}{1};
            if ~isempty(EEG) & ~isempty(tY_table{r,c}{1})
                tY = tY_table{r,c}{1}; Y = tY(:,:,2); t = tY(:,:,1);
                evs = EEG.event; 
                evs = evs(arrayfun(@(ev) ~sum(strcmp(ev.type, {'-1','15','14','12','13'})), evs));
                init_time = [evs.init_time];
                timeDiff = diff(init_time);
                sameTrial = [0, timeDiff <= timeBetweenEvents];
                % ^ event i = same trial as event i-1
                startEvs = evs(~sameTrial);
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
                    for evIdx = 1:length(curEEG.event)
                        curEEG.event(evIdx).latency = curEEG.event(evIdx).latency + startT*curEEG.srate;
                    end
                    EEGs(idx) = curEEG;

                    ti = (t >= startT) & (t < startT+maxTrialDur);
                    curY = Y; curY(~ti) = nan;
                    curT = t; curT(~ti) = nan;
                    tYs{idx} = cat(3,curT,curY);

                    clear startT curEEG ti curY curT evIdx
                end
                EEG_trial{r,c} = {EEGs};
                tY_trial{r,c} = {tYs};
            end
            clear EEG tY init_time timeDiff sameTrial EEGs tYs evs startEvs
        end
    end
    trialTables{subj,1} = tY_trial; trialTables{subj,2} = EEG_trial;
    clear EEG_table tY_table EEG_trial tY_trial
end

%% organize for plotting 
for v = testVars
    disp(['setting up plots for ',v{:}]);

    disp('counting all trials')
    trialCount = zeros(size(trialTables,1),1);
    for subj = 1:size(trialTables,1)
        EEG_trial = trialTables{subj,2}; EEG_trial = makeSubtbl(EEG_trial, v);
        tY_trial  = trialTables{subj,1}; tY_trial  = makeSubtbl(tY_trial,  v);
        for r = 1:height(EEG_trial)
            trialCount(subj) = trialCount(subj) + ...
                min(length(EEG_trial{r,1}{:}), length(tY_trial{r,1}{:}));
        end
    end
    H = sum(trialCount > 0); %trialTables = trialTables((trialCount > 0), :);
    W = max(trialCount);
    W = W+1; % plotting channel colors 

    disp('plotting time series')
    figure; sgtitle([yname,' ',v{:}]);
    idx3 = 0; 
    for subj = 1:length(scanfiles)
        idx3incr = false;
        fn = scanfiles{subj}
        strend = find(fn == '-'); strend = strend((diff(strend)==1)); strend = max(1, strend(1)-2);
        ylbl = fn(1:strend);
        ylbl = ylbl((end-2):end);

        BL = BLs{subj,1};
        EEG_trial = trialTables{subj,2}; EEG_trial = makeSubtbl(EEG_trial, v);
        tY_trial  = trialTables{subj,1}; tY_trial  = makeSubtbl(tY_trial,  v);
        idx1 = 1; idx2 = 1;
        for r = 1:height(EEG_trial)
            EEGs = EEG_trial{r,1}{:};
            tYs = tY_trial{r,1}{:};
            for trl = 1:length(EEGs)
                if ~isempty(EEGs) & ~isempty(tYs)
                    EEG = EEGs(trl);
                    tY = tYs{trl};
                    if ~isempty(EEG) & ~isempty(tY)
                        chlocs = EEG.chanlocs;
                        Y = tY(:,:,2); t = tY(:,:,1);
                        chanSig = sum( ( (Y < BL(1,:)) | (Y > BL(3,:)) ),1 );

                        idx2 = (idx3)*W + idx1; idx3incr = true;
                        subplot(H,W,idx2);
                        ttl = [EEG_trial.Properties.RowNames{r},' trial ',num2str(trl)];
                        plotEvents(EEG, ylims);
                        title(ttl); ylim(ylims); ylabel(ylbl); xlabel('time (s)');

                        for chan = 1:size(Y,2)
                            if chanSig(chan)
                                tc = t(:,chan); Yc = Y(:,chan);
                                plot(tc, Yc, ...
                                    'Color',chanColor(chlocs(chan),chlocs),...
                                    'LineWidth',1);
                            end
                        end
                        %legend({chlocs.labels})
                        %{
                        for chan = 1:size(Y,2)
                            if chanSig(chan)
                                tc = t(:,chan); Yc = Y(:,chan); 
                                fill([min(tc),min(tc),max(tc),max(tc)],...
                                     [BL(1,chan),BL(3,chan),BL(3,chan),BL(1,chan)],...
                                     chanColor(chlocs(chan), chlocs),...
                                     'FaceAlpha',baselineTransparency);
                            end
                        end
                        %}

                        idx1 = idx1 + 1;
                    end
                end
            end
        end
        if idx3incr
            subplot(H,W,idx2+1); hold on;
            for chan = chlocs
                plot3(chan.X, chan.Y, chan.Z, '.', ...
                    'Color', chanColor(chan, chlocs));
                text(chan.X, chan.Y, chan.Z, chan.labels, ...
                    'Color', chanColor(chan, chlocs));
            end
            %idx1 = idx1 + 1;
        end
        idx3 = idx3 + idx3incr;
    end
    clear idx1 idx2 EEG_trial tY_trial EEGs tYs EEG tY t Y ttl ylbl strend ...
          chanSig chan chlocs tc Yc idx3 idx3incr;

    disp('running hypothesis tests')
    W = W-1;
    figure; sgtitle([yname,' ',v{:}]);
    for subj = 1:H
        fn = scanfiles{subj}
        strend = find(fn == '-'); strend = strend((diff(strend)==1)); strend = max(1, strend(1)-2);
        ylbl = fn(1:strend);
        ylbl = ylbl((end-2):end);

        EEG_trial = trialTables{subj,2}; EEG_trial = makeSubtbl(EEG_trial, v);
        tY_trial  = trialTables{subj,1}; tY_trial  = makeSubtbl(tY_trial,  v);
        BL = BLs{subj,2}; BL_t = BL(:,:,1); BL = BL(:,:,2);
        idx1 = 1;
        for r = 1:height(EEG_trial)
            EEGs = EEG_trial{r,1}{:};
            tYs = tY_trial{r,1}{:};
            for trl = 1:length(EEGs)
                if ~isempty(EEGs) & ~isempty(tYs)
                    EEG = EEGs(trl);
                    tY = tYs{trl};
                    if ~isempty(EEG) & ~isempty(tY)
                        Y = tY(:,:,2); t = tY(:,:,1);

                        pp = zeros(1,size(Y,2));
                        for chan = 1:size(Y,2)
                            y = Y(:,chan); bl = BL(:,chan);
                            y = y(~isnan(y)); bl = bl(~isnan(bl));
                            if mean(y) > mean(bl)
                                [~,pp(chan)] = ttest2(bl,y, 'Vartype','unequal', 'Tail','left');
                                pp(chan) = 1-pp(chan);
                            else
                                [~,pp(chan)] = ttest2(bl,y, 'Vartype','unequal', 'Tail','right');
                                pp(chan) = pp(chan)-1;
                            end
                        end

                        idx2 = (subj-1)*W + idx1;
                        subplot(H,W,idx2);
                        topoplot(pp, EEG.chanlocs, 'maplimits', [-1,1]); colorbar;
                        title([EEG_trial.Properties.RowNames{r},' trial ',num2str(trl)]);
                        ylabel(ylbl);

                        idx1 = idx1 + 1;
                    end
                end
            end
        end
    end 
    clear idx1 indx2 EEG_trial tY_trial EEGs tYs EEG tY Y_t Y BL BL_t...
          pp ttl ylbl strend chan y bl;
end

%% helper functions 

function subtbl = makeSubtbl(tbl, vars)
    subtbl = tbl(:, ismember(tbl.Properties.VariableNames, vars));
end

function rgb = chanColor(chloc, chlocs)
    if isempty(chloc.X)
        chloc.X = 0;
    end
    if isempty(chloc.Y)
        chloc.Y = 0;
    end
    if isempty(chloc.Z)
        chloc.Z = 0;
    end
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
                eegTbl{r,c} = {eegTbl{r,c}{1}(1)};
            end
        end
    end
end


function plotEvents(EEG, ybound)
    event = EEG.event; srate = EEG.srate;
    hold on; 
    lbl_sw = false;
    for ev = event
        if ~isempty(ev.init_time) & ~strcmp(ev.type,'boundary')
            initTime = ev.latency/srate;
            %initTime = ev.init_time;
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