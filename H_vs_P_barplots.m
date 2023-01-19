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
svloc = [postproDir,'/Poster Stats ',...
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

figure; hold on;
[~,idx] = unique(upper({allchan.labels})); 
allchan = allchan(idx);
for chan = allchan
    plot3(chan.X, chan.Y, chan.Z, '.', ...
        'Color', chanColor(chan, allchan));
    text(chan.X, chan.Y, chan.Z, chan.labels, ...
        'Color', chanColor(chan, allchan));
end
[chansel, chanselName] = listdlg_selectWrapper({allchan.labels}, ...
    'multiple', 'Select Channels:');
allchan0 = allchan; allchan = allchan(chansel);

%% combination "subjects" 
comboSubj = cell(2, length(scanfiles));
DATATABLES_chansel = DATATABLES;

for s = 1:length(scanfiles)

    somechan = true(size(allchan));
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

    somechan = allchan(somechan);
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
                tY = tY(:,ord,:);
                cumuTY = cat(1, cumuTY, tY);
                tY_trial{trl} = tY;
            end
            cumuTYsubj{c} = cumuTY;
            dataTable{4,c} = {tY_trial};
        end
        cumuTYs = [cumuTYs; cumuTYsubj];
        dataTables{subj} = dataTable;
    end

    cumuTYsubj = cell(1, size(cumuTYs,2));
    for c = 1:size(cumuTYs,2)
        cumuTY = [];
        for subj = 1:size(cumuTYs,1)
            cumuTY = cat(1, cumuTY, cumuTYs{subj, c});
        end
        cumuTYsubj{c} = cumuTY;
    end

    comboSubj{1, s} = cumuTYsubj; comboSubj{2, s} = somechan;
    DATATABLES_chansel{s} = dataTables;
    clear cumuTY cumuTYs cumuTYsubj somechan dataTable dataTables ...
          EEG_all tY_trial tY ord 
end

%% plots
maxplt = -inf(size(chansel)); minplt = inf(size(chansel));

%clr = {'b', 'r', 'k', 'm', 'c', 'g', 'y'};
clr = {[  1, 66,130]/255, ... navy blue 
       [127,105,  3]/255, ... darkened yellow 
       [ 17,172,228]/255, ... light blue
       [254,209,  6]/255, ... yellow
       [247,132, 34]/255, ... orange 
       [ 66,174, 73]/255, ... green 
       [213, 32, 39]/255, ... red
       [141,153,193]/255, ... purple blue 
       [244,154,192]/255, ... pink
       };
mkr = {'o', '^', 's', 'p', 'h', 'd', '>', 'v', '<'};
spc2 = 2; spc1 = 3; spc3 = 4;

brplt = zeros(size(comboSubj,2), ...
        max( length(comboSubj{1,1}),length(comboSubj{1,2}) ), ...
        length(chanselName), ...
        2);

fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(yname);
for s = 1:length(scanfiles)
    sf = scanfiles{s};
    dataTables = DATATABLES{s};

    for subj = 1:length(sf)
        fn = sf{subj}; pname = fn(1:3);
        dataTable = dataTables{subj};

        for c = 1:width(dataTable)
            v = dataTable.Properties.VariableNames{c};
            tY_trial = dataTable{4, c}{1}; % extra wrapping shouldn't be necessary?
            EEG_all  = dataTable{1, c}{1};

            Y_val = zeros(length(tY_trial), size(tY_trial{1},2));
            Y_erb = Y_val;
            for trl = 1:length(tY_trial)
                tY = tY_trial{trl};
                Y = tY(:,:,2);
                Y_val(trl,:) = mean(Y, 'omitnan');
                Y_erb(trl,:) = std(Y, 'omitnan'); % change to SE or 95% CI?
            end

            xplt = length(scanfiles)*(maxNgrp*(maxNtrl+spc3)+spc1)*spc2*(c-1) + ...
                   (maxNgrp*(maxNtrl+spc3)+spc1)*(s-1) + ...
                   (subj-1)*(maxNtrl+spc3);

            for ch = 1:length(chanselName)
                subplot(length(chanselName),1,ch); grid on; hold on;

                chName = chanselName{ch};
                chIdx = strcmpi(chName, {EEG_all.chanlocs.labels});

                yplt = Y_val(:,chIdx); yplte = Y_erb(:,chIdx);

                ylabel(chName);
                xtk = (length(scanfiles)*(maxNgrp*(maxNtrl+spc3)+spc1)*spc2*...
                    ((1:width(dataTable))-1) + maxNgrp*(maxNtrl+spc3));
                xticks(xtk);
                xticklabels(dataTable.Properties.VariableNames);
                xlim([-spc2*(maxNgrp*(maxNtrl+spc3)+spc1), ...
                    length(scanfiles)*(maxNgrp*(maxNtrl+spc3)+spc1)*spc2*width(dataTable)]);

                errorbar(xplt+(1:length(yplt))-1, ...
                    yplt, yplte, ...
                    'Color',clr{s}, 'Marker',mkr{subj}, 'LineStyle','none', 'LineWidth',.5);
            end

        end
    end

    % cumulative "subject" 
    cumuTYs = comboSubj{1, s}; cumuchan = comboSubj{2, s};
    for c = 1:width(cumuTYs)
        tY = cumuTYs{c}; 
        Y = tY(:,:,2);
        Y_val = mean(Y, 'omitnan');
        Y_erb = std(Y, 'omitnan'); % change to SE or 95% CI?

        xplt = length(scanfiles)*(maxNgrp*(maxNtrl+spc3)+spc1)*spc2*(c-1) + ...
               (maxNgrp*(maxNtrl+spc3)+spc1)*(s-1) + ...
               (subj)*(maxNtrl+spc3);

        for ch = 1:length(chanselName)
            subplot(length(chanselName),1,ch); 

            chName = chanselName{ch};
            chIdx = strcmpi(chName, {cumuchan.labels});

            yplt = Y_val(:,chIdx); yplte = Y_erb(:,chIdx);
            maxplt(ch) = max(maxplt(ch), max(yplt + yplte));
            minplt(ch) = min(minplt(ch), min(yplt - yplte));

            errorbar(xplt, yplt, yplte, ...
                'Color',clr{s}, 'Marker','*', 'LineStyle','none', 'LineWidth',2);

            brplt(s,c,ch,1) = yplt; brplt(s,c,ch,2) = yplte;
        end

    end

end
clear sf dataTables v fn pname tY_trial EEG_all ...
    Y_val Y_erb tY Y yplt yplte xplt chName chIdx ...
    cumuchan cumuTYs

for ch = 1:length(chanselName)
    subplot(length(chanselName),1,ch);
    ax = gca; ax.FontSize = 16;
end

saveas(fig, [svloc,yname,' AllTrialsPlot'], 'fig');

%% export baseline stat as mat file 
MeasurementName = yname;
BaselineMeasurementTable = table('RowNames', chanselName);

for s = 1:length(scanfiles)
    sf = scanfiles{s};
    dataTables = DATATABLES_chansel{s};
    for subj = 1:length(sf)
        fn = sf{subj} 
        pname = fn(1:3);
        dataTable = dataTables{subj};

        blTable = dataTable.BaselineBefore;
        tY_trial = blTable{4};
        EEG_all = blTable{1};
        if length(tY_trial) > 1
            warning('Using first baseline trial only.')
        end
        Ystat = struct(...
            'mean', 0, ...
            'SD',   0, ...
            'SE',   0);
            Y_stats = repmat(Ystat, length(tY_trial), size(tY_trial{1},2));
            trl = 1;
            tY = tY_trial{trl}; Y = tY(:,:,2);
            for ch = 1:size(tY,2)
                y = Y(:,ch);
                Y_stats(ch).mean = mean(y, 'omitnan');
                Y_stats(ch).SD   = std(y, 'omitnan'); 
                Y_stats(ch).SE   = Y_stats(ch).SD / sqrt(length(y));
            end
            Y_stats = Y_stats';
        eval(['BaselineMeasurementTable.',pname,' = Y_stats']);
    end
end

%% hypothesis testing 
% make more robust with more than 2 experimental groups / more baselines? 

pvals = zeros(3, ...
        max( length(comboSubj{1,1}),length(comboSubj{1,2}) ), ...
        length(chanselName) );
tYs = comboSubj(1,:);
for s = 1:length(tYs)
    cumuchan = comboSubj{2,s};
    cumuchansel = false(size(cumuchan));
    for ch = 1:length(cumuchansel)
        cumuchansel(ch) = sum(strcmpi(cumuchan(ch).labels, chanselName));
    end
    tYs_s = tYs{s};
    for c = 1:length(tYs_s)
        tYs_s{c} = tYs_s{c}(:,cumuchansel,:);
    end
    tYs{s} = tYs_s;
end
for c = 1:size(pvals,2)
    for s = 1:2
        tYs_s = tYs{s};
        for ch = 1:size(pvals,3)
            [~,pvals(s,c,ch)] = ...
                ttest2( tYs_s{1}(:,ch,2), ...
                        tYs_s{c}(:,ch,2), ...
                        'Vartype', 'unequal' );
        end
    end
    for ch = 1:size(pvals,3)
        [~,pvals(3,c,ch)] = ...
            ttest2( tYs{1}{c}(:,ch,2), ...
                    tYs{2}{c}(:,ch,2), ...
                    'Vartype', 'unequal' );
    end
end
clear tYs tYs_s cumuchan cumuchansel

fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(yname);
ytk = 3; alph = 0.001;
spch = .15; 
for ch = 1:length(chansel)
    subplot(length(chanselName),1,ch);
    avg = brplt(:,:,ch,1)'; 
    b = bar(avg); hold on; grid on;
    for s = 1:length(b)
        b(s).FaceColor = clr{s+2};
    end
    for s = 1:2
        errorbar((1:size(brplt,2))+(s-1.5)*2*spch, ...
            brplt(s,:,ch,1), brplt(s,:,ch,2), ...
            'Color','k', 'LineStyle','none', 'LineWidth',2);
    end
    xticklabels(dataTable.Properties.VariableNames);
    ylabel(chanselName{ch});
    spcv = (maxplt(ch) - minplt(ch))/ytk;
    % below
    curY = minplt(ch) - spcv;
    for c = 1:size(pvals,2)
        if pvals(3,c,ch) < alph
            plotPval(pvals(3,c,ch), c-spch, c+spch, curY);
            %curY = curY - spcv;
        end
    end
    % above 
    curY = maxplt(ch) + spcv;
    for c = 1:size(pvals,2)
        if pvals(1,c,ch) < alph
            plotPval(pvals(1,c,ch), 1-spch, c-spch, curY);
            curY = curY + spcv;
        end
        if pvals(2,c,ch) < alph
            plotPval(pvals(2,c,ch), 1+spch, c+spch, curY);
            curY = curY + spcv;
        end
    end
    ylim([minplt(ch)-2*spcv, curY]);
end
for ch = 1:length(chanselName)
    subplot(length(chanselName),1,ch);
    ax = gca; ax.FontSize = 16;
end

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

function plotPval(p, x1, x2, y)
    x = mean([x1, x2]); xbar = abs(x2-x1)./2;
    errorbar(x,y,0,0,xbar,xbar, 'k', 'LineWidth',1);
    text(x,y, ['p = ',num2str(p,1)], 'FontSize',8, 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom');
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