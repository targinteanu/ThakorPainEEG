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
scanfiles = scanfiles((~strcmp({scanfiles.name},'.')) & ...
                      (~strcmp({scanfiles.name},'..')) & ...
                      (~strcmp({scanfiles.name},'.DS_Store')) );
scanfiles = scanfiles(~[scanfiles.isdir]);
[~,scanfilesH] = listdlg_selectWrapper({scanfiles.name},'multiple','Select Controls');
[~,scanfilesP] = listdlg_selectWrapper({scanfiles.name},'multiple','Select Patients');
scanfiles = {scanfilesH, scanfilesP};
scanfileNames = {'Control', 'Patient'};
maxNgrp = max( arrayfun(@(s) length(scanfiles{s}), 1:length(scanfiles)) );
maxNgrp = maxNgrp + 1; % make room for combo subj

cd(home); addpath(postproDir);
svloc = [postproDir,'/Comparison Topoplots ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
mkdir(svloc); svloc = [svloc,'/'];

%% select what to plot

[fcn, yname, ylims] = MeasurementSelector();

%% calculations 
timeBetweenEvents = 5; timeAfterLast = 3; % seconds 
DATATABLES = cell(size(scanfiles));

for s = 1:length(scanfiles)
    sf = scanfiles{s};

    dataTables = cell(length(sf),1);
    for subj = 1:length(sf)
        fn = sf{subj}
        load(fn);

        % run calculations on desired variables
        disp('calculating desired variables')
        % construct data tables *************
        tempTbl = table('RowNames',{'EEG_all','EEG_trial','tY_all','tY_trial'});

        tempArr = EEG_table.BaselineOpen('before experiment');
        tempTbl.BaselineBefore('EEG_all') = tempArr; 
        curSpec = EpochSpec_table.BaselineOpen('before experiment'); curSpec = curSpec{:};
        curEpoc = Epoch_table.BaselineOpen('before experiment');     curEpoc = curEpoc{:};
        [Y,t] = fcn(curSpec, curEpoc); 
        tempTbl.BaselineBefore('tY_all') = {cat(3,t,Y)}; 

        tempArr = EEG_table.BaselineOpen('after experiment');
        tempTbl.BaselineAfter('EEG_all') = tempArr; 
        curSpec = EpochSpec_table.BaselineOpen('after experiment'); curSpec = curSpec{:};
        curEpoc = Epoch_table.BaselineOpen('after experiment');     curEpoc = curEpoc{:};
        [Y,t] = fcn(curSpec, curEpoc); 
        tempTbl.BaselineAfter('tY_all') = {cat(3,t,Y)}; 

        tempArr1 = EEG_table.TempStim('before experiment'); tempArr2 = EEG_table.TempStim('after experiment');
        tempArr = [tempArr1{:} tempArr2{:}];
        tempTbl.TempStim('EEG_all') = {tempArr}; 
        tempArr1 = Epoch_table.TempStim('before experiment'); tempArr2 = Epoch_table.TempStim('after experiment');
        curEpoc = ArrCat(tempArr1{:}, tempArr2{:});
        tempArr1 = EpochSpec_table.TempStim('before experiment'); tempArr2 = EpochSpec_table.TempStim('after experiment');
        curSpec = ArrCat(tempArr1{:}, tempArr2{:});
        [Y,t] = fcn(curSpec, curEpoc);
        tempTbl.TempStim('tY_all') = {cat(3,t,Y)}; 

        tempArr1 = EEG_table.PinPrick('before experiment'); tempArr2 = EEG_table.PinPrick('after experiment');
        tempArr = [tempArr1{:} tempArr2{:}];
        tempTbl.PinPrick('EEG_all') = {tempArr}; 
        tempArr1 = Epoch_table.PinPrick('before experiment'); tempArr2 = Epoch_table.PinPrick('after experiment');
        curEpoc = ArrCat(tempArr1{:}, tempArr2{:});
        tempArr1 = EpochSpec_table.PinPrick('before experiment'); tempArr2 = EpochSpec_table.PinPrick('after experiment');
        curSpec = ArrCat(tempArr1{:}, tempArr2{:});
        [Y,t] = fcn(curSpec, curEpoc);
        tempTbl.PinPrick('tY_all') = {cat(3,t,Y)}; 

        %{
        [tY_table, EEG_table] = ...
            fcnTbl(EEG_table, Epoch_table, EpochSpec_table, fcn, testVars);
        %}

        dataTables{subj} = tempTbl;
        clear tempTbl tempArr tempArr1 tempArr2 curEpoc curSpec
        clear tempTbl EEG_table Epoch_table EpochSpec_table Spec_table
    end

    DATATABLES{s} = dataTables;
    clear sf dataTables
end

%% determine max trial duration
disp('determining trial durations')
maxTrialDur = 0; % s
for DT = 1:length(DATATABLES)
    dataTables = DATATABLES{DT};
    for subj = 1:length(dataTables)
        EEG_table = dataTables{subj}{1,:};
        for c = 1:length(EEG_table)
            EEG = EEG_table{c};
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
        clear EEG_table
    end
    clear dataTables
end

%% segment into trials 
maxNtrl = 0;
for DT = 1:length(DATATABLES)
    dataTables = DATATABLES{DT};
    sf = scanfiles{DT};
    for subj = 1:length(dataTables)
        fn = sf{subj}
        dataTable = dataTables{subj};

        EEG_table = dataTable{1,:};
        EEG_trial = EEG_table;
        tY_table  = dataTable{3,:};
        tY_trial  = tY_table;

        for c = 1:length(EEG_table)
            v = dataTable.Properties.VariableNames{c};
            disp(['segmenting ',v,' into trials'])

            for r = 1:height(EEG_table)
                EEG = EEG_table{c};
                if ~isempty(EEG) & ~isempty(tY_table{c})
                    % setup variables
                    tY = tY_table{c}; Y = tY(:,:,2); t = tY(:,:,1);
                    evs = EEG.event;
                    % don't consider events of unwanted type
                    evs = evs(arrayfun(@(ev) ~sum(strcmp(ev.type, {'-1','15','14','12','13'})), evs));

                    % --------------------------------------------------------

                    if sum(strcmp(v, {'BaselineBefore', 'BaselineAfter'}))
                        startEvs = evs(1); endEvs = evs(end);
                        startEvs = [startEvs.latency]/EEG.srate;
                        endEvs = [endEvs.latency]/EEG.srate;

                    elseif strcmp(v, 'PinPrick')
                        evs = evs(~~arrayfun(@(ev) sum(strcmp(ev.type, {'11','10'})), evs));

                        init_time = [evs.init_time];
                        timeDiff = diff(init_time);
                        intvlFromPrev = [inf, timeDiff]; intvlToNext = [timeDiff, inf];

                        t_inter_ev = 2; t_before_ev = 10;
                        firstOfTrain = (intvlToNext < t_inter_ev) & (intvlFromPrev >= t_inter_ev);
                        lastOfTrain = (intvlFromPrev < t_inter_ev) & (intvlToNext >= t_inter_ev);
                        prickBefore = (intvlToNext >= t_inter_ev) & (intvlToNext < t_before_ev) & (intvlFromPrev >= t_inter_ev);

                        evStart = find(firstOfTrain); evEnd = zeros(size(evStart));
                        for idx = 1:length(evStart)
                            eEnd = find(lastOfTrain);
                            eEnd = eEnd(eEnd >= evStart(idx));
                            if idx < length(evStart)
                                eEnd = eEnd(eEnd <= evStart(idx+1));
                            end
                            if isempty(eEnd)
                                evEnd(idx) = -1;
                            else
                                evEnd(idx) = max(eEnd);
                            end
                        end
                        evStart = evStart(evEnd >= 0); evEnd = evEnd(evEnd >= 0);

                        ev0 = zeros(size(evStart));
                        for idx = 1:length(evStart)
                            e0 = find(prickBefore);
                            e0 = e0(e0 <= evStart(idx));
                            if idx > 1
                                e0 = e0(e0 >= evStart(idx-1));
                            end
                            if isempty(e0)
                                ev0(idx) = -1;
                            else
                                ev0(idx) = min(e0);
                            end
                        end

                        startEvs = ev0;
                        startEvs(ev0 < 0) = evStart(ev0 < 0);
                        startEvs = evs(startEvs); endEvs = evs(evEnd);
                        startEvs = [startEvs.latency]/EEG.srate;
                        endEvs = [endEvs.latency]/EEG.srate;

                    elseif strcmp(v, 'TempStim')
                        evs = evs(arrayfun(@(ev) strcmp(ev.type, '3'), evs));

                        init_time = [evs.init_time];
                        timeDiff = diff(init_time);
                        intvlFromPrev = [inf, timeDiff]; intvlToNext = [timeDiff, inf];

                        t_inter_ev = 4;
                        firstOfTrain = (intvlToNext < t_inter_ev) & (intvlFromPrev >= t_inter_ev);
                        lastOfTrain = (intvlFromPrev < t_inter_ev) & (intvlToNext >= t_inter_ev);

                        evStart = find(firstOfTrain); evEnd = zeros(size(evStart));
                        for idx = 1:length(evStart)
                            eEnd = find(lastOfTrain);
                            eEnd = eEnd(eEnd >= evStart(idx));
                            if idx < length(evStart)
                                eEnd = eEnd(eEnd <= evStart(idx+1));
                            end
                            if isempty(eEnd)
                                evEnd(idx) = -1;
                            else
                                evEnd(idx) = max(eEnd);
                            end
                        end
                        evStart = evStart(evEnd >= 0); evEnd = evEnd(evEnd >= 0);

                        startEvs = evs(evStart); endEvs = evs(evEnd);
                        startEvs = [startEvs.latency]/EEG.srate;
                        endEvs = [endEvs.latency]/EEG.srate;

                    elseif strcmp(v, 'Pressure')
                        startEvs = evs(~~arrayfun(@(ev) sum(strcmp(ev.type, {'4','8'})), evs));
                        endEvs   = evs(~~arrayfun(@(ev) sum(strcmp(ev.type, {'5','9'})), evs));

                        % order by time
                        if length(startEvs) > 1
                            [~,ord] = sort([startEvs.init_time]); startEvs = startEvs(ord);
                        end
                        % match start/end
                        pressEv = [];
                        for startEv_ = startEvs
                            endEvOpts = endEvs([endEvs.init_time] >= startEv_.init_time);
                            if ~isempty(endEvOpts)
                                [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                                endEv_ = endEvOpts(endEv_);
                                pressEv = [pressEv; startEv_, endEv_];
                            end
                        end
                        if isempty(pressEv)
                            startEvs = []; endEvs = [];
                        else
                            startEvs = pressEv(:,1); endEvs = pressEv(:,2);

                            startEvs = [startEvs.latency]/EEG.srate;
                            endEvs = [endEvs.latency]/EEG.srate;
                        end

                    else

                        init_time = [evs.init_time];
                        timeDiff = diff(init_time);
                        intvlFromPrev = [inf, timeDiff]; intvlToNext = [timeDiff, inf];

                        sameTrial = [0, timeDiff <= timeBetweenEvents];
                        % ^ event i = same trial as event i-1
                        startEvs = evs(~sameTrial);
                        startEvs = [startEvs.latency]/EEG.srate;
                        endEvs = startEvs + maxTrialDur;

                    end

                    endEvs = endEvs + timeAfterLast;

                    clear t_inter_ev t_before_ev firstOfTrain lastOfTrain prickBefore
                    clear intvlToNext intvlFromPrev timeDiff init_time
                    clear pressEv startEv_ endEv_ endEvOpts ord

                    EEGs = repmat(EEG, size(startEvs));
                    tYs = cell(size(startEvs));

                    for idx = 1:length(startEvs)
                        startT = startEvs(idx);
                        endT = endEvs(idx);

                        disp(['Segment ',num2str(idx),' of ',num2str(length(startEvs)),...
                            ' (',num2str(100*idx/length(startEvs),3),'%)'])
                        curEEG = extractBetweenTimes(EEG, [startT, endT]);
                        curEEG.xmin = curEEG.xmin + startT;
                        curEEG.xmax = curEEG.xmax + startT;
                        curEEG.times = curEEG.times + startT*1000;
                        for evIdx = 1:length(curEEG.event)
                            curEEG.event(evIdx).latency = curEEG.event(evIdx).latency + startT*curEEG.srate;
                        end
                        EEGs(idx) = curEEG;

                        ti = (t >= startT) & (t < endT);
                        ti = ti | [zeros(1,size(ti,2));ti(1:(end-1),:)] | [ti(2:end,:);zeros(1,size(ti,2))] ;
                        curY = Y; curY(~ti) = nan;
                        curT = t; curT(~ti) = nan;
                        tYs{idx} = cat(3,curT,curY);

                        clear startT endT curEEG ti curY curT evIdx
                    end
                    maxNtrl = max(maxNtrl, length(startEvs));
                    EEG_trial{r,c} = EEGs;
                    tY_trial{r,c} = tYs;
                    % --------------------------------------------------------

                end
                clear EEG tY init_time timeDiff sameTrial EEGs tYs evs startEvs
            end
        end
        
        dataTable{2,:} = EEG_trial; dataTable{4,:} = tY_trial;
        dataTables{subj} = dataTable; 
        clear EEG_table tY_table EEG_trial tY_trial
    end
    DATATABLES{DT} = dataTables;
    clear dataTables dataTable sf
end

%% channel selection 
%%{
% get all channel names of all EEGs 
allchan = [];
for s = 1:length(scanfiles)
    sf = scanfiles{s};
    dataTables = DATATABLES{s};
    for subj = 1:length(sf)
        dataTable = dataTables{subj};
        for c = 1:width(dataTable)
            EEG = dataTable{1,c}{1};
            allchan = [allchan, EEG.chanlocs];
        end
    end
end
clear sf dataTable dataTables EEG

% eliminate repeats (case insensitive) 
[~,idx] = unique(upper({allchan.labels})); 
allchan = allchan(idx);

% only include chans common to all EEGs 
idx = true(size(allchan));
for s = 1:length(scanfiles)
    sf = scanfiles{s};
    dataTables = DATATABLES{s};
    for subj = 1:length(sf)
        dataTable = dataTables{subj};
        for c = 1:width(dataTable)
            EEG = dataTable{1,c}{1};
            for ch = 1:length(allchan)
                idx(ch) = idx(ch) & ...
                    sum( strcmpi(allchan(ch).labels, {EEG.chanlocs.labels}) );
            end
        end
    end
end
allchan = allchan(idx);
clear sf dataTable dataTables EEG idx

% select desired chans 
figure; hold on;
for chan = allchan
    plot3(chan.X, chan.Y, chan.Z, '.', ...
        'Color', chanColor(chan, allchan));
    text(chan.X, chan.Y, chan.Z, chan.labels, ...
        'Color', chanColor(chan, allchan));
end
[chansel, chanselName] = listdlg_selectWrapper({allchan.labels}, ...
    'multiple', 'Select Channels:');
allchan0 = allchan; allchan = allchan(chansel);
%}
%% combination "subjects" 
varnames = DATATABLES{1}{1}.Properties.VariableNames;
comboSubjTbl = table('size',[length(scanfileNames),length(varnames)], ...
                     'VariableTypes',repmat("cell",size(varnames)), ...
                     'VariableNames',varnames,'RowNames',scanfileNames);

somechan = true(size(allchan));
for s = 1:length(scanfiles)
    dataTables = DATATABLES{s};
    for subj = 1:length(dataTables)
        dataTable = dataTables{subj};
        for c = 1:width(dataTable)
            EEG_all = dataTable{1,c}{1};
            for ch = 1:length(somechan)
                somechan(ch) = somechan(ch) & sum( ...
                    strcmpi(allchan(ch).labels, {EEG_all.chanlocs.labels}) );
            end
        end
    end
end
somechan = allchan(somechan);
% somechan = only those selected AND in all EEGs of all subjects of all groups

for s = 1:length(scanfiles)

    dataTables = DATATABLES{s};
    cumuTYs = {};
    for subj = 1:length(dataTables)
        dataTable = dataTables{subj};
        cumuTYsubj = cell(1, width(dataTable));
        for c = 1:width(dataTable)
            EEG_all  = dataTable{1,c}{1};
            tY_trial = dataTable{4,c}{1};
            cumuTY = [];
            for trl = 1:length(tY_trial)
                tY = tY_trial{trl};
                ord = zeros(1,length(somechan));
                for ch = 1:length(somechan)
                    ord(ch) = find( ...
                        strcmpi(somechan(ch).labels, {EEG_all.chanlocs.labels}) );
                end
                cumuTY = cat(1, cumuTY, tY(:,ord,:)); % all trials concatenated 
            end
            cumuTYsubj{c} = cumuTY; % diff stimuli within one subj
        end
        cumuTYs = [cumuTYs; cumuTYsubj]; % #subjs x #stimtypes
    end

    % concat all subjs (collapse cumuTYs vertically)
    for c = 1:size(cumuTYs,2)
        cumuTY = [];
        for subj = 1:size(cumuTYs,1)
            cumuTY = cat(1, cumuTY, cumuTYs{subj, c});
        end
        comboSubjTbl{s,c} = {cumuTY};
    end

    clear cumuTY cumuTYs cumuTYsubj dataTable dataTables ...
          EEG_all tY_trial tY ord 
end
% everything in comboSubj should be ordered according to somechan

%% comparison 
% make more robust with more than 2 experimental groups / more baselines? 
varnames = DATATABLES{1}{1}.Properties.VariableNames;

statsTable = table('size',[3,length(varnames)], ...
                   'VariableTypes',repmat("cell",size(varnames)), ...
                   'RowNames',{'P vs H','H vs baseline','P vs baseline'}, ...
                   'VariableNames',varnames);

p_alpha = 0.01; % uncertainty for stat significance 
maxstatval = -Inf; minstatval = Inf;

% statistical testing 
for c = 1:width(comboSubjTbl)

    s = repmat(struct('chan',somechan(1), 'tstat',0, 'pval',1), ...
               [length(somechan),height(statsTable)]);
    for ch = 1:length(somechan)
        % P vs H
        [~,p,~,S] = ...
            ttest2( comboSubjTbl{2,c}{1}(:,ch,2), ...
                    comboSubjTbl{1,c}{1}(:,ch,2), ...
                    'Vartype', 'unequal' );
        s(ch,1).chan  = somechan(ch);
        s(ch,1).tstat = S.tstat;
        s(ch,1).pval  = p;

        % H vs baseline 
        [~,p,~,S] = ...
            ttest2( comboSubjTbl{1,c}{1}(:,ch,2), ...
                    comboSubjTbl{1,1}{1}(:,ch,2), ...
                    'Vartype', 'unequal' );
        s(ch,2).chan  = somechan(ch);
        s(ch,2).tstat = S.tstat;
        s(ch,2).pval  = p;

        % P vs baseline 
        [~,p,~,S] = ...
            ttest2( comboSubjTbl{2,c}{1}(:,ch,2), ...
                    comboSubjTbl{2,1}{1}(:,ch,2), ...
                    'Vartype', 'unequal' );
        s(ch,3).chan  = somechan(ch);
        s(ch,3).tstat = S.tstat;
        s(ch,3).pval  = p;

        clear p S
    end

    statsTable{1,c} = {s(:,1)};
    statsTable{2,c} = {s(:,2)};
    statsTable{3,c} = {s(:,3)};
    clear s
end

% for plotting: find bounds and significance 
for r = 1:height(statsTable)
    for c = 1:width(statsTable)
        S = statsTable{r,c}{1};
        %{
        for ch = 1:length(S)
            if S(ch).pval < p_alpha
                S(ch).chan.labels = '*';
            else
                S(ch).chan.labels = '.';
            end
        end
        %}
        statsTable{r,c} = {S};
        S = [S.tstat];
        maxstatval = max(maxstatval, max(S)); minstatval = min(minstatval, min(S));
    end
end
clear S
maxstatval = max(abs(maxstatval), abs(minstatval));
minstatval = -maxstatval;

%% plot comparison 
plot_vs_baseline = {'TempStim', 'PinPrick'};
plot_P_vs_H      = {'BaselineBefore'};

fig = figure('Units', 'Normalized', 'Position', [0 0 1 .3]); 
sgtitle([yname,' t statistic; * p < ',num2str(p_alpha)]);
W = 2*length(plot_vs_baseline) + length(plot_P_vs_H);

w = 1;
for idx = 1:length(plot_P_vs_H)
    subplot(1,W, w);
    pltname = plot_P_vs_H{idx};
    S = statsTable(1, strcmp(pltname, statsTable.Properties.VariableNames));
    title([pltname,' ',S.Properties.RowNames{1}]);
    S = S{1,1}{1};
    topoplot([S.tstat], [S.chan], ...
        'emarker2', {find([S.pval]<p_alpha), '*', 'k'}, ...
        ...'pmask', [S.pval]<p_alpha, ...
        ...'numcontour',[p_alpha p_alpha], 'contourvals',[S.pval], ...
        'maplimits', [minstatval, maxstatval]); colorbar;
    w = w + 1;
end
for idx = 1:length(plot_vs_baseline)
    pltname = plot_vs_baseline{idx};
    S = statsTable([2,3], strcmp(pltname, statsTable.Properties.VariableNames));
    for idx2 = 1:height(S)
        subplot(1,W, w);
        title([pltname,' ',S.Properties.RowNames{idx2}]);
        SS = S{idx2,1}{1};
        topoplot([SS.tstat], [SS.chan], ...
            'emarker2', {find([SS.pval]<p_alpha), '*', 'k'}, ...
            'maplimits', [minstatval, maxstatval]); colorbar;
        w = w + 1;
    end
end
clear S SS

saveas(fig, [svloc,yname,' StatsBarPlot'], 'fig');

%% helper functions 

function outArr = ArrCat(Arr1, Arr2)
    if isempty(Arr1)
        outArr = Arr2;
    elseif isempty(Arr2)
        outArr = Arr1;
    else
        outArr = [Arr1, Arr2];
    end
end

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

function outEEG = extractBetweenTimes(inEEG, bound)
    buf = .5; % s
    bound(1) = max(inEEG.xmin, bound(1)-buf);
    bound(2) = min(inEEG.xmax, bound(2)+buf);
    outEEG = pop_select(inEEG, 'time', bound);
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