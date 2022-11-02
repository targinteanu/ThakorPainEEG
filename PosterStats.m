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

cd(home); addpath(postproDir);
svloc = [postproDir,'/Poster Stats ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
mkdir(svloc); svloc = [svloc,'/'];

%% select what to plot

[fcn, yname, ylims] = MeasurementSelector();

[~,testVars] = listdlg_selectWrapper(...
    {'TempStim', 'PinPrick', 'Pressure', 'BaselineOpen', 'BaselineClosed', 'BaselineIce'}, ...
    'multiple', 'Specify Variable');

%% calculations 
timeBetweenEvents = 5; timeAfterLast = 3; % seconds 
DATATABLES = cell(size(scanfiles));

for s = 1:length(scanfiles)
    sf = scanfiles{s};

    dataTables = cell(length(sf),2);
    for subj = 1:length(sf)
        fn = sf{subj}
        load(fn);

        % run calculations on desired variables
        disp('calculating desired variables')
        [tY_table, EEG_table] = ...
            fcnTbl(EEG_table, Epoch_table, EpochSpec_table, fcn, testVars);

        dataTables{subj,1} = tY_table; dataTables{subj,2} = EEG_table;
        clear tY_table EEG_table Epoch_table EpochSpec_table Spec_table
    end

    DATATABLES{s} = dataTables;
    clear sf dataTables
end

%% determine max trial duration 
disp('determining trial durations')
maxTrialDur = 0; % s
for DT = 1:length(DATATABLES)
    dataTables = DATATABLES{DT};
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
    clear dataTables
end

%% segment into trials 
TRIALTABLES = cell(size(DATATABLES));
for DT = 1:length(DATATABLES)
    dataTables = DATATABLES{DT};
    sf = scanfiles{DT};
    trialTables = cell(size(dataTables));
    for subj = 1:size(dataTables,1)
        fn = sf{subj}

        EEG_table = dataTables{subj,2};
        EEG_trial = EEG_table;
        tY_table  = dataTables{subj,1};
        tY_trial  = tY_table;

        for v = testVars
            v = v{:};
            disp(['segmenting ',v,' into trials'])

            c = find(strcmp(v, EEG_table.Properties.VariableNames));

            for r = 1:height(EEG_table)
                EEG = EEG_table{r,c}{1};
                if ~isempty(EEG) & ~isempty(tY_table{r,c}{1})
                    % setup variables
                    tY = tY_table{r,c}{1}; Y = tY(:,:,2); t = tY(:,:,1);
                    evs = EEG.event;
                    % don't consider events of unwanted type
                    evs = evs(arrayfun(@(ev) ~sum(strcmp(ev.type, {'-1','15','14','12','13'})), evs));

                    % --------------------------------------------------------

                    if strcmp(v, 'BaselineOpen')
                        startEvs = evs(1); endEvs = evs(end);
                        startEvs = [startEvs.latency]/EEG.srate;
                        endEvs = [endEvs.latency]/EEG.srate;

                    elseif strcmp(v, 'BaselineClosed')
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
                    EEG_trial{r,c} = {EEGs};
                    tY_trial{r,c} = {tYs};
                    % --------------------------------------------------------

                end
                clear EEG tY init_time timeDiff sameTrial EEGs tYs evs startEvs
            end
        end
        trialTables{subj,1} = tY_trial; trialTables{subj,2} = EEG_trial;
        clear EEG_table tY_table EEG_trial tY_trial
    end
    TRIALTABLES{DT} = trialTables;
    clear dataTables sf trialTables
end


%% plots 
clr = {'b', 'r', 'k', 'm', 'c', 'g', 'y'};
mkr = {'o', 's', 'p', 'h', 'x', 'd', '^', '>', 'v', '<', '+'};
spc = 10*maxNgrp; 

for s = 1:length(scanfiles)
    sf = scanfiles{s};
    dataTables = DATATABLES{s};
    trialTables = TRIALTABLES{s};

    for subj = 1:length(sf)
        fn = sf{subj};
        tY_table = trialTables{subj, 1};

        for c = 1:width(tY_table)
            for r = 1:height(tY_table)
                tY_trial = tY_table{r, c}{1};
                v = tY_table.Properties.VariableNames{c};
                cond = tY_table.Properties.RowNames{r};

                for trl = 1:length(tY_trial)
                    tY = tY_trial{trl};
                    Y = tY(:,:,2);
                    Y_val = mean(Y, 'omitnan');
                    Y_erb = std(Y, 'omitnan'); % change to SE or 95% CI?
                end

            end
        end
    end
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