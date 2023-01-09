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

%% table reorganization - necessary? 
timeBetweenEvents = 5; timeAfterLast = 3; % seconds 
DATATABLES = cell(size(scanfiles));

for s = 1:length(scanfiles)
    sf = scanfiles{s};

    dataTables = cell(length(sf),1);
    for subj = 1:length(sf)
        fn = sf{subj}
        load(fn);

        % construct data tables *************
        tempTbl = table('RowNames',{'EEG_all','EEG_trial','Epocs','Specs','Spec_trial','Spec_alltrial'});
        tempArr = EEG_table.BaselineOpen('before experiment');
        tempTbl.BaselineBefore('EEG_all') = tempArr; 
        curSpec = EpochSpec_table.BaselineOpen('before experiment'); %curSpec = curSpec{:};
        curEpoc = Epoch_table.BaselineOpen('before experiment');     %curEpoc = curEpoc{:};
        tempTbl.BaselineBefore('Epocs') = curEpoc; tempTbl.BaselineBefore('Specs') = curSpec;

        tempArr = EEG_table.BaselineOpen('after experiment');
        tempTbl.BaselineAfter('EEG_all') = tempArr; 
        curSpec = EpochSpec_table.BaselineOpen('after experiment'); %curSpec = curSpec{:};
        curEpoc = Epoch_table.BaselineOpen('after experiment');     %curEpoc = curEpoc{:};
        tempTbl.BaselineAfter('Epocs') = curEpoc; tempTbl.BaselineAfter('Specs') = curSpec;

        tempArr1 = EEG_table.TempStim('before experiment'); tempArr2 = EEG_table.TempStim('after experiment');
        tempArr = [tempArr1{:} tempArr2{:}];
        tempTbl.TempStim('EEG_all') = {tempArr}; 
        tempArr1 = Epoch_table.TempStim('before experiment'); tempArr2 = Epoch_table.TempStim('after experiment');
        curEpoc = ArrCat(tempArr1{:}, tempArr2{:});
        tempArr1 = EpochSpec_table.TempStim('before experiment'); tempArr2 = EpochSpec_table.TempStim('after experiment');
        curSpec = ArrCat(tempArr1{:}, tempArr2{:});
        tempTbl.TempStim('Epocs') = curEpoc; tempTbl.TempStim('Specs') = curSpec;

        tempArr1 = EEG_table.PinPrick('before experiment'); tempArr2 = EEG_table.PinPrick('after experiment');
        tempArr = [tempArr1{:} tempArr2{:}];
        tempTbl.PinPrick('EEG_all') = {tempArr}; 
        tempArr1 = Epoch_table.PinPrick('before experiment'); tempArr2 = Epoch_table.PinPrick('after experiment');
        curEpoc = ArrCat(tempArr1{:}, tempArr2{:});
        tempArr1 = EpochSpec_table.PinPrick('before experiment'); tempArr2 = EpochSpec_table.PinPrick('after experiment');
        curSpec = ArrCat(tempArr1{:}, tempArr2{:});
        tempTbl.PinPrick('Epocs') = curEpoc; tempTbl.PinPrick('Specs') = curSpec; 

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
        Epocs  = dataTable{3,:};
        Specs  = dataTable{4,:};
        Specs_trial = Specs;

        for c = 1:length(EEG_table)
            v = dataTable.Properties.VariableNames{c};
            disp(['segmenting ',v,' into trials'])

            for r = 1:height(EEG_table)
                EEG = EEG_table{c};
                Spec = Specs{c}; Epoc = Epocs{c};
                if ~isempty(EEG) & ~isempty(Spec) & ~isempty(Epoc)
                    % setup variables
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
                    Spec_tf = false(size(Spec)); 
                    Spec_trial = cell(size(startEvs));

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

                        for idx2 = 1:length(Spec)
                            curEpoc = Epoc(idx2);
                            Spec_tf(idx2) = ...
                                min(curEpoc.times) <=   endT*1000    & ... change to | for more permissive inclusion
                                max(curEpoc.times) >= startT*1000;
                        end
                        Spec_trial{idx} = Spec(Spec_tf);

                        clear startT endT curEEG ti evIdx idx2 Spec_tf curEpoc
                    end
                    maxNtrl = max(maxNtrl, length(startEvs));
                    EEG_trial{r,c} = EEGs;
                    Specs_trial{r,c} = Spec_trial;
                    % --------------------------------------------------------

                end
                clear EEG init_time timeDiff sameTrial EEGs evs startEvs Spec_trial Spec Epoc
            end
        end
        
        dataTable{2,:} = EEG_trial; dataTable{5,:} = Specs_trial;
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


%% reorganize channels for concatenation and ignore unselected 
varnames = DATATABLES{1}{1}.Properties.VariableNames;

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

w_all = [];
for s = 1:length(scanfiles)
    dataTables = DATATABLES{s};
    for subj = 1:length(dataTables)
        dataTable = dataTables{subj};
        for c = 1:width(dataTable)
            EEG_all = dataTable{1,c}{1};
            ord = zeros(1,length(somechan));
            for ch = 1:length(somechan)
                ord(ch) = find( ...
                    strcmpi(somechan(ch).labels, {EEG_all.chanlocs.labels}) );
            end
            Specs = dataTable{5,c}{1};
            Specs_alltrial = [];
            for trl = 1:length(Specs)
                Spec = Specs{trl};
                for idx = 1:length(Spec)
                    Spec_i = Spec(idx);
                    Spec_i.frequencySpectrum = Spec_i.frequencySpectrum(ord,:);
                    Spec_i.powerSpectrum     = Spec_i.powerSpectrum(ord,:);
                    Spec_i.chanlocs          = Spec_i.chanlocs(ord);

                    w_all = [w_all, Spec_i.frequency1side];
                    w_all = unique(w_all);

                    Spec(idx) = Spec_i;
                    clear Spec_i
                end

                % collapse trials for each subject 
                Specs_alltrial = [Specs_alltrial, Spec];

                Specs{trl} = Spec; 
                clear Spec
            end
            dataTable{5,c} = {Specs}; dataTable{6,c} = {Specs_alltrial};
            clear Specs Specs_alltrial EEG_all ord
        end
        dataTables{subj} = dataTable; 
        clear dataTable
    end
    DATATABLES{s} = dataTables; 
    clear dataTables
end
w_all = sort(unique(w_all));

%% collapse subjects for each subject group 
comboSubjTbl = table('size',[length(scanfileNames),length(varnames)], ...
                     'VariableTypes',repmat("cell",size(varnames)), ...
                     'VariableNames',varnames,'RowNames',scanfileNames);

for s = 1:length(scanfiles)
    dataTables = DATATABLES{s};

    for c = 1:width(varnames)
        P_all = [];
        for subj = 1:length(dataTables)
            dataTable = dataTables{subj};

            EEG_all = dataTable{1,c}{1};
            Specs = dataTable{5,c}{1};
            Specs_alltrial = dataTable{6,c}{1};

            P2 = zeros(length(somechan), length(w_all), length(Specs_alltrial));
            for idx = 1:length(Specs_alltrial)
                Spec_i = Specs_alltrial(idx);
                P = Spec_i.powerSpectrum; w = Spec_i.frequency1side;
                P2(:,:,idx) = interp1(w', P', w_all', 'linear', 'extrap')';
                clear Spec_i P w 
            end
            P_all = cat(3, P_all, P2);
        end
        clear dataTable

        comboSpecs.frequency1side = w_all;
        comboSpecs.powerSpectrum = P_all;
        comboSpecs.chanlocs = somechan;
        comboSpecs.nbchan = length(somechan);
        comboSubjTbl{s,c} = {comboSpecs};

    end
    clear dataTables
end

%% plotting 


%% helper functions 

function outArr = ArrCat(Arr1, Arr2)
    if isempty(Arr1)
        outArr = Arr2;
    elseif isempty(Arr2)
        outArr = Arr1;
    else
        outArr = [Arr1, Arr2];
    end
    outArr = {outArr};
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