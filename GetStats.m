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
mkdir(svloc); svloc = [svloc,'/'];

%% select what to plot
meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);
plotOpts = {'Peak Freq (Hz)', 'Band Mean (\muV^2 s^2)', 'Band Density', ...
            'Node Degree', 'Connectivity Strength', 'Frequency Assortativity', ...
            'Neighbors Average '};
analysisOpts = {'single variable', 'correlation'};
analysisType = questdlg('Calculate What?', 'Analysis Type', ...
                        analysisOpts{1}, analysisOpts{2}, analysisOpts{1});
analysisType = find(strcmp(analysisType, analysisOpts));
if analysisType == 1
    nvar = 'single';
    yname = [];
elseif analysisType == 2
    nvar = 'multiple';
    yname = 'Correlation Between ';
end
fcnOpts = {@(w,P,band) peakFreq(w,P,band), meanAmp, ampDensity};
PLOTSEL = listdlg_selectWrapper(plotOpts, nvar, 'Plot What?');
if ~(length(PLOTSEL) == analysisType)
    if (analysisType == 2) & (length(PLOTSEL) == 1)
        PLOTSEL = [PLOTSEL, PLOTSEL];
    else
        error('Incorrect number of selections.')
    end
end

fcn = cell(1,2);
for idx = 1:length(PLOTSEL)
plotSel = PLOTSEL(idx);
neighborLayers = 0;
if sum(plotSel == [4,6,7])
    % percentile cutoff must be selected 
    Pcut = inputdlg('Cutoff Percentile (%)','Network Cutoff Selection');
    Pcut = str2num(Pcut{1});
end
while plotSel == length(plotOpts)
    % Neighbors Average Value - of what?
    neighborLayers = neighborLayers + 1;
    yname = [yname, plotOpts{plotSel}];
    plotSel = listdlg_selectWrapper(plotOpts, 'single', 'Average What?');
end
if sum(plotSel == 1:6)
    % frequency band must be selected 
    [bnd,bndname] = pickFrequency();
    yname = [yname,bndname,' ',plotOpts{plotSel}];
    if sum(plotSel == [2,3])
        fcn0 = fcnOpts{plotSel};
        bnd = band2freqs(bnd, BandTableHz);
        fcn{idx} = @(Spec, EEG) frqFcnEpoch(Spec, EEG, @(w,P) fcn0(w,P,bnd));
        if plotSel == 3
            % Density 
            ylims = [0 1];
        else
            ylims = [];
        end
    elseif plotSel == 1
        % peak freq
        fcn{idx} = @(Spec, EEG) peakFreqEpoch(Spec, EEG, bnd, BandTableHz);
        if isa(bnd,'char') | isa(bnd,'string')
            ylims = band2freqs(bnd, BandTableHz);
        else
            ylims = bnd;
        end
    elseif plotSel == 6
        % assortativity 
        fcn{idx} = @(Spec, EEG) assortativity(Spec, EEG, bnd, Pcut, BandTableHz);
        ylims = [-1, 1];
    elseif plotSel == 4
        % node degree
        fcn{idx} = @(Spec, EEG) nodeDegree(Spec, EEG, Pcut, bnd, BandTableHz);
        ylims = 'numchan';
    elseif plotSel == 5
        % conn strength
        fcn{idx} = @(Spec, EEG) connStrength(Spec, EEG, bnd, BandTableHz);
        ylims = [];
    end
end
for l = 1:neighborLayers
    fcn{idx} = @(Spec, EEG) avg_neighbor(Spec, EEG, fcn{idx}, Pcut, bnd, BandTableHz);
end

if analysisType == 2
    if idx == 1
        yname = [yname,' and '];
    else
        fcn = @(Spec, EEG) fcnCorr(fcn{1}, fcn{2}, Spec, EEG);
        ylims = [-1, 1];
    end
else
    fcn = fcn{1};
end
end

[~,testVars] = listdlg_selectWrapper(...
    {'TempStim', 'PinPrick', 'Pressure', 'BaselineOpen', 'BaselineClosed', 'BaselineIce'}, ...
    'multiple', 'Specify Variable');

%% baseline and calculations  
%table2baseline = @(tbl) [tbl.BaselineOpen('before experiment'),...
%                         tbl.BaselineOpen('after experiment')];
table2baseline = @(tbl) tbl.BaselineOpen('before experiment');
p = .05; % uncertainty level 
timeBetweenEvents = 5; timeAfterLast = 15; % seconds 
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
            [Y,t] = fcn(BL_EpochSpec{idx}, BL_Epoch{idx});
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
    title(['Subject ',pname,' Baseline']);
    if ~isempty(ylims)
        if strcmp(ylims, 'numchan')
            templims = [0,size(BL,2)];
        else
            templims = ylims;
        end
        ylim(templims);
    else
        templims = 'maxmin';
    end

    if length(BL(2,:)) > 1
        xticklabels({EEG.chanlocs.labels}); xlabel('Channel'); ylabel(yname);
        
        figure(fig(2)); sgtitle(['Baseline ',yname]);
        subplot(W,H,subj); 
        topoplot(BL(2,:), EEG.chanlocs, 'maplimits', templims, 'electrodes','labels'); colorbar;
        title(['Subject ',pname]);
    end
end
saveas(fig(1), [svloc,yname,'_Baseline_Bar_Plot'], 'fig');
saveas(fig(2), [svloc,yname,'_Baseline_Head_Map'], 'fig');
clear fig BL EEG EEGrc W H templims

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
    fig = figure; sgtitle([yname,' ',v{:}]);
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
                        
                        if ~isempty(ylims)
                            if strcmp(ylims, 'numchan')
                                templims = [0,size(BL,2)];
                            else
                                templims = ylims;
                            end
                            ylim(templims);
                        else
                            templims = [min([BL(:);Y(:)]), max([BL(:);Y(:)])];
                        end

                        plotEvents(EEG, templims);
                        title(ttl); ylim(templims); ylabel(ylbl); xlabel('time (s)');

                        for chan = 1:size(Y,2)
                            if chanSig(chan)
                                tc = t(:,chan); Yc = Y(:,chan);
                                plot(tc, Yc, ...
                                    'Color',chanColor(chlocs(chan),chlocs),...
                                    'LineWidth',1);
                            end
                        end
                        %legend({chlocs.labels})
                        %%{
                        for chan = 1:size(Y,2)
                            if chanSig(chan)
                                tc = t(:,chan); Yc = Y(:,chan); 
                                fill([min(tc),min(tc),max(tc),max(tc)],...
                                     [BL(1,chan),BL(3,chan),BL(3,chan),BL(1,chan)],...
                                     chanColor(chlocs(chan), chlocs),...
                                     'FaceAlpha',baselineTransparency, 'EdgeColor','none');
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
    saveas(fig, [svloc,yname,'_',v{:},'_Trial_Time_Series'], 'fig');
    clear idx1 idx2 EEG_trial tY_trial EEGs tYs EEG tY t Y ttl ylbl strend ...
          chanSig chan chlocs tc Yc idx3 idx3incr;

    BL = BLs{1,2};  
    if size(BL,2) > 1
    disp('running hypothesis tests')
    W = W-1;
    fig = figure; sgtitle([yname,' ',v{:}]);
    idx3 = 0;
    for subj = 1:length(scanfiles)
        idx3incr = false;
        fn = scanfiles{subj}
        strend = find(fn == '-'); strend = strend((diff(strend)==1)); strend = max(1, strend(1)-2);
        ylbl = fn(1:strend);
        ylbl = ylbl((end-2):end);

        EEG_trial = trialTables{subj,2}; EEG_trial = makeSubtbl(EEG_trial, v);
        tY_trial  = trialTables{subj,1}; tY_trial  = makeSubtbl(tY_trial,  v);
        BL = BLs{subj,2}; BL_t = BL(:,:,1); BL = BL(:,:,2);
        idx1 = 1; idx2 = 1;
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
                                %[~,pp(chan)] = ttest2(bl,y, 'Vartype','unequal', 'Tail','left');
                                [~,pp(chan)] = kstest2(bl,y, 'Tail','larger');
                                pp(chan) = .5-pp(chan);
                            else
                                %[~,pp(chan)] = ttest2(bl,y, 'Vartype','unequal', 'Tail','right');
                                [~,pp(chan)] = kstest2(bl,y, 'Tail','smaller');
                                pp(chan) = pp(chan)-.5;
                            end
                        end

                        idx2 = (idx3)*W + idx1; idx3incr = true;
                        subplot(H,W,idx2);
                        topoplot(pp, EEG.chanlocs, 'maplimits', [-.5,.5]); colorbar;
                        title([EEG_trial.Properties.RowNames{r},' trial ',num2str(trl)]);
                        ylabel(ylbl);

                        idx1 = idx1 + 1;
                    end
                end
            end
        end
        idx3 = idx3 + idx3incr;
    end 
    saveas(fig, [svloc,yname,'_',v{:},'_Trial_Head_Map'], 'fig');
    end
    clear idx1 indx2 EEG_trial tY_trial EEGs tYs EEG tY Y_t Y BL BL_t...
          pp ttl ylbl strend chan y bl;
end

%% helper functions 

function subtbl = makeSubtbl(tbl, vars)
    subtbl = tbl(:, ismember(tbl.Properties.VariableNames, vars));
end

function t = getTimes(epoch_EEG)
    if isempty(epoch_EEG)
        t = [];
    else
        t = arrayfun(@(eeg) mean([eeg.xmin, eeg.xmax]), epoch_EEG);
    end
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

function [Y,t] = fcnCorr(fcn1, fcn2, var1, var2)
    [Y1,t1] = fcn1(var1, var2); [Y2,t2] = fcn2(var1, var2);
    Y1 = Y1'; Y2 = Y2'; t1 = t1'; t2 = t2';

    t = sort(unique([t1(:);t2(:)]));
    Y1 = cell2mat( arrayfun(@(c) ...
        interp1(t1(c,:), Y1(c,:), t, 'linear', 'extrap'), ...
        1:size(Y1,1), 'UniformOutput',false) )';
    Y2 = cell2mat( arrayfun(@(c) ...
        interp1(t2(c,:), Y2(c,:), t, 'linear', 'extrap'), ...
        1:size(Y2,1), 'UniformOutput',false) )';

    Y = zeros(length(t),1);
    for s = 1:length(t)
        Y(s) = corr(Y1(:,s), Y2(:,s), 'Type', 'Spearman', 'Rows', 'complete');
    end
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
                [Y,t] = fcn(curSpec, curEpoc);
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

%% time-frequency functions 

function [Y,times] = frqFcnEpoch(epoch_Spec, epoch_EEG, fcn)
    Y = cell2mat( arrayfun(@(eegSpec) ...
        fcn(eegSpec.frequency1side, eegSpec.powerSpectrum), ...
        epoch_Spec, 'UniformOutput', false) );
    times = getTimes(epoch_EEG);
    times = repmat(times,size(Y,1),1);
    times = times'; Y = Y';
end

function [Freq,times] = peakFreqEpoch(epoch_Spec, epoch_EEG, bnd, tbl)
% inspect peakFreq across epochs (time)
    nchan = epoch_Spec(1).nbchan;
    PAFs = zeros(nchan, length(epoch_Spec), 2);
    for chan = 1:nchan
        PAFs(chan,:,1) = arrayfun(@(eegSpec) ...
            peakFreq(eegSpec.frequency1side, eegSpec.powerSpectrum(chan,:), ...
            bnd, tbl), ...
            epoch_Spec);
        PAFs(chan,:,2) = arrayfun(@(eeg) ...
            mean([eeg.xmin, eeg.xmax]), epoch_EEG);
    end
    times = PAFs(:,:,2)'; Freq = PAFs(:,:,1)';
end

%% network functions 

function [nd,t] = nodeDegree(SpectObj, EEGObj, cutoffPercentile, bnd, tbl)
    if nargin < 5
        tbl = [];
        if nargin < 4
            bnd = [];
            if nargin < 3
                cutoffPercentile = [];
            end
        end
    end
    nd = zeros(SpectObj(1).nbchan,length(SpectObj));
    for s = 1:length(SpectObj)
        [W,A] = WPLI(SpectObj(s).frequencySpectrum, cutoffPercentile, ...
            bnd, SpectObj(s).frequency2side, tbl);
        nd(:,s) = sum(A);
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(nd,1), 1); 
    nd = nd'; t = t';
end

function [cs,t] = connStrength(SpectObj, EEGObj, bnd, tbl)
    if nargin < 4
        tbl = [];
        if nargin < 3
            bnd = [];
        end
    end
    % each channel's average WPLI with all others 
    cs = zeros(SpectObj(1).nbchan,length(SpectObj));
    for s = 1:length(SpectObj)
        W = WPLI(SpectObj(s).frequencySpectrum, [], bnd, SpectObj(s).frequency2side, tbl);
        cs(:,s) = mean(W, 'omitnan');
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(cs,1), 1); 
    cs = cs'; t = t';
end

function [A,t] = avg_neighbor(SpectObj, EEGObj, neighborFcn, cutoffPercentile, bnd, tbl)
    if nargin < 6
        tbl = [];
        if nargin < 5
            bnd = [];
            if nargin < 4
                cutoffPercentile = [];
            end
        end
    end
    Y = neighborFcn(SpectObj, EEGObj); Y = Y';
    A = zeros(size(Y));
    for s = 1:length(SpectObj)
        [~,Adj] = WPLI(SpectObj(s).frequencySpectrum, cutoffPercentile, bnd, ...
            SpectObj(s).frequency2side, tbl);
        A(:,s) = arrayfun(@(c) mean(Y(Adj(:,c),s)), 1:size(Y,1));
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(A,1), 1); 
    A = A'; t = t';
end

function [Af,t] = assortativity(SpectObj, EEGObj, bnd, cutoffPercentile, tbl)
    % 1D
    if nargin < 5
        if ~exist('BandTableHz')
            load('BrainwaveFrequencyTable.mat');
            global BandTableHz
        end
        tbl = BandTableHz;
        if nargin < 4
            cutoffPercentile = [];
            if nargin < 3
                bnd = [];
                if nargin < 2
                    EEGObj = [];
                end
            end
        end
    end
    Af = zeros(size(SpectObj));
    for s = 1:length(SpectObj)
        [~,Adj] = WPLI(SpectObj(s).frequencySpectrum, cutoffPercentile, bnd, ...
            SpectObj(s).frequency2side, tbl);
        [~,PF] = diffFreq(SpectObj(s).frequency1side, SpectObj(s).powerSpectrum, bnd, tbl);
        y = arrayfun(@(c) mean(PF(Adj(:,c))), 1:length(PF));
        Af(s) = corr(PF', y', 'Type','Spearman', 'Rows', 'complete');
    end
    t = getTimes(EEGObj);
    Af = Af'; t = t';
end