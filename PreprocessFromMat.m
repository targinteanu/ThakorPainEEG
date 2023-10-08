%% Set Filepaths
clear

% functions for manipuliating dirs:
% get only subfolders
subfoldersof = @(d) d([d.isdir] & ...
    ~strcmp({d.name}, '.') & ...
    ~strcmp({d.name}, '..'));

%home = 'C:\Users\targi\Desktop\Thakor Chronic Pain Data\Data_Chronic Pain';
home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)
datafolders = dir;
datafolders = subfoldersof(datafolders);

% user selects which scans to use
[sel, ok] = listdlg('ListString',{datafolders.name}, 'SelectionMode','multiple');
while ~ok
    sel = questdlg('select all?');
    ok = strcmp(sel, 'Yes');
    if ~ok
        [sel, ok] = listdlg('ListString',{datafolders.name}, 'SelectionMode','multiple');
    else
        sel = 1:length(datafolders);
    end
end
datafolders = datafolders(sel);

%addpath 'C:\Users\targi\Documents\MATLAB\ThakorPainEEG';
addpath = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
svloc = [home,'/Preprocessed ',datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];

%% Start eeglab
eeglabpath = 'C:\Program Files\MATLAB\R2022a\eeglab2023.0';
addpath(eeglabpath)
eeglab

%% Preprocessing in EEGLAB ---------------------------------------------
mkdir(svloc);
%%
for subj = 1:size(datafolders,1)
    cd(home)
    disp(datafolders(subj,1).name)
    cd(datafolders(subj,1).name);

        clearvars datasets eventsets
        datasets = dir('*.mat');
        N = size(datasets,1);

        for d2 = 1:N
            clearvars EEG rejected EEG_table ...
                pcaComp pcaScr pcaPexp rejectedPCA ...
                endEv endEv_ startEv startEv_ endEvOpts eventType eventTime ord toExclude ...
                iceEv iceEvInit iceEvFin openEv openEvInit openEvFin closedEv closedEvInit closedEvFin ...
                pressEv pressEvInit pressEvFin pressCpmEv prickEv prickEvInit prickEvFin prickCpmEv ...
                tempEv tempEvInit tempEvFin 

            id = [datafolders(subj,1).name,' --- ',datasets(d2,1).name];

            % import EEG and EV2
            load(datasets(d2,1).name);
            EEG = eeg_checkset( EEG );
            EEG = pop_chanedit(EEG, 'lookup', [eeglabpath,'/plugins/dipfit/standard_BEM/elec/standard_1005.elc']);
            EEG = eeg_checkset( EEG );

            %%{

            % initial modification 
            if EEG.srate > 512
                EEG = pop_resample( EEG, 512);
            end
            EEG = pop_select( EEG,'nochannel',{'CB1' 'CB2' 'HEO' 'VEO' 'HEOG' 'VEOG' 'M1' 'M2'});
            EEG = pop_reref( EEG, []);
            EEG = pop_eegfiltnew(EEG, 1, 100, 826, 0, [], 1);
            notch = designfilt('bandstopiir', ...
                'PassbandFrequency1', 58, 'StopbandFrequency1', 59, ...
                'StopbandFrequency2', 61, 'PassbandFrequency2', 62, ...
                'PassbandRipple1', 1, 'PassbandRipple2', 1, 'StopbandAttenuation', 10, ...
                'SampleRate', EEG.srate);
            EEG.data = single(filtfilt(notch, double(EEG.data'))');
            EEG = eeg_checkset(EEG);

            %inspect raw spectra and reject spectra that are odd looking
            figure; pop_spectopo(EEG, 1, [0  15000], 'EEG' , 'freqrange',[2 100],'electrodes','off');
            pop_eegplot( EEG, 1, 1, 1);
            rejected = input('Please identify sensor # to remove: ');
            if size(rejected,1)
                EEG = pop_select( EEG,'nochannel',rejected);
                EEG = pop_reref( EEG, []); % KD - reref again after removing noisy channels
            end

            % remove eye movement artifacts 
            [pcaComp, pcaScr, ~,~, pcaPexp] = pca(EEG.data');
            nrows = floor(sqrt(size(pcaComp,2))); ncols = ceil(size(pcaComp,2)/nrows);
            figure; 
            for idx = 1:size(pcaComp,2)
                subplot(nrows, ncols, idx);
                topoplot(pcaComp(:,idx), EEG.chanlocs);
                title(['PC ',num2str(idx)])
            end
            eegplot(pcaScr');
            rejectedPCA = input('Please identify PCA component # to remove: ');
            cont = ~isempty(rejectedPCA); cont0 = cont;
            while cont
                rejector = ones(size(pcaScr,2),1); rejector(rejectedPCA) = 0;
                rejector = diag(rejector);
                outData = (pcaScr*rejector)*(pcaComp^-1); outData = outData';
                eegplot(outData); 
                sel = questdlg('Proceed with removed components?', ...
                    'PCA confirmation', 'Yes', 'Retry', 'Cancel', 'Retry');
                cont = strcmp(sel, 'Retry') | strcmp(sel, '');
                if cont
                    eegplot(pcaScr');
                    rejectedPCA = input('Please identify PCA component # to remove: ');
                    cont = ~isempty(rejectedPCA);
                end
            end
            if cont0
                if strcmp(sel, 'Yes')
                    EEG.data = outData;
                    EEG = pop_reref( EEG, []); 
                end
            end
            clearvars cont cont0 outData nrows ncols idx



            %}


            % ================ EXTRACT DESIRED EVENTS ====================
            eventType = [EEG.event.type]; eventTime = [EEG.event.init_time];

            % BASELINE ---------------------------------------------------
            
            % eyes open: 0 to 1
                startEv = EEG.event(eventType == 0); endEv = EEG.event(eventType == 1);
                excludeTypes = 1:13; % exclude baseline if any of these is found between start/end
                splitTypes = 3:13; % END if any of these is found before baseline start
                [openEv, openEvInit, openEvFin] = getBaselineEvs(startEv, endEv, excludeTypes, splitTypes);

            % baseline eyes closed: 1 to 2
                startEv = EEG.event(eventType == 1); endEv = EEG.event(eventType == 2);
                excludeTypes = [0, 2:13]; % exclude baseline if any of these is found between start/end 
                splitTypes = 3:13; % END if any of these is found before baseline start
                [closedEv, closedEvInit, closedEvFin] = getBaselineEvs(startEv, endEv, excludeTypes, splitTypes);

            % hand in ice: 6 to 7
                startEv = EEG.event(eventType == 6); endEv = EEG.event(eventType == 7);
                excludeTypes = [6, 7]; % exclude baseline if any of these is found between start/end
                iceEv = getBaselineEvs(startEv, endEv, excludeTypes, []);

            % TEMP STIM --------------------------------------------------

            % temp stim: 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                startEv = EEG.event(eventType == 3); 
                excludeTypes = [-1, 3, 5, 9, 12:15]; % these types DO NOT end interval !!!
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                tempEv = matchStartEndEvs(startEv, endEv);
                % designate before/after ice experiment 
                splitTypes = [6, 7, 8:10, 12, 13];
                [tempEvInit, tempEvFin] = splitBeforeAfter(tempEv, splitTypes);

            % PIN PRICK --------------------------------------------------

            % pin prick: 11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                startEv = EEG.event(eventType == 11); 
                excludeTypes = [-1, 11, 5, 9, 12:15]; % these types DO NOT end interval !!!
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                prickEv = matchStartEndEvs(startEv, endEv);
                % designate before/after ice experiment 
                splitTypes = [6, 7, 8:10, 12, 13];
                [prickEvInit, prickEvFin] = splitBeforeAfter(prickEv, splitTypes);

            % CPM pin prick: 10 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                startEv = EEG.event(eventType == 10); 
                excludeTypes = [-1, 10, 5, 9, 12:15]; % these types DO NOT end interval !!!
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                prickCpmEv = matchStartEndEvs(startEv, endEv);
                % select only those within iceEv????

            % PRESSURE ---------------------------------------------------

            % pressure: 4 (to 5)  
                startEv = EEG.event(eventType == 4); 
                endEv = EEG.event(eventType == 5); 
                pressEv = matchStartEndEvs(startEv, endEv);
                % designate before/after ice experiment 
                splitTypes = [6, 7, 8:10, 12, 13];
                [pressEvInit, pressEvFin] = splitBeforeAfter(pressEv, splitTypes);

            % CPM pressure: 8 (to 9)  
                startEv = EEG.event(eventType == 8); 
                endEv = EEG.event(eventType == 9); 
                pressCpmEv = matchStartEndEvs(startEv, endEv);
                % select only those within iceEv????

            % ------------------------------------------------------------

            % baseline in ice
            allEv = [prickCpmEv; pressCpmEv; prickEv; pressEv; tempEv];
            if isempty(allEv)
                iceEvInit = iceEv; iceEvFin = [];
            else
                iceEvInit = []; iceEvFin = [];
                for idx = 1:size(iceEv,1)
                    evOpts = ([allEv(:,2).init_time] > iceEv(idx,1).init_time) & ...
                             ([allEv(:,1).init_time] < iceEv(idx,2).init_time);
                    evOpts = allEv(evOpts,:);
                    [ev1time,ev1idx] = min([evOpts(:,1).init_time]);
                    if ev1time > iceEv(idx,1).init_time
                        iceEvInit = [iceEvInit; iceEv(idx,1), evOpts(ev1idx,1)];
                    end
                    for startEv_ = evOpts(:,2)'
                        endEvOpts = evOpts(([evOpts(:,2).init_time] > startEv_.init_time),:);
                        [ev1time,ev1idx] = min([endEvOpts(:,1).init_time]);
                        if ev1time > startEv_.init_time
                            iceEvFin = [iceEvFin; startEv_, endEvOpts(ev1idx,1)];
                        end
                    end
                end
            end

            % ============================================================

            % extract EEG from events 
            EEG_table = table('RowNames',{'before experiment','after experiment','CPM'});
            EEG_table.BaselineOpen(1) = extractBetweenEvents(EEG, uniqueEvents(openEvInit));
            EEG_table.BaselineOpen(2) = extractBetweenEvents(EEG, uniqueEvents(openEvFin));
            EEG_table.BaselineClosed(1) = extractBetweenEvents(EEG, uniqueEvents(closedEvInit));
            EEG_table.BaselineClosed(2) = extractBetweenEvents(EEG, uniqueEvents(closedEvFin));
            EEG_table.BaselineIce(1) = extractBetweenEvents(EEG, uniqueEvents(iceEvInit));
            EEG_table.BaselineIce(2) = extractBetweenEvents(EEG, uniqueEvents(iceEvFin));
            EEG_table.TempStim(1) = extractBetweenEvents(EEG, uniqueEvents(tempEvInit));
            EEG_table.TempStim(2) = extractBetweenEvents(EEG, uniqueEvents(tempEvFin));
            EEG_table.PinPrick(1) = extractBetweenEvents(EEG, uniqueEvents(prickEvInit));
            EEG_table.PinPrick(2) = extractBetweenEvents(EEG, uniqueEvents(prickEvFin));
            EEG_table.PinPrick(3) = extractBetweenEvents(EEG, uniqueEvents(prickCpmEv));
            EEG_table.Pressure(1) = extractBetweenEvents(EEG, uniqueEvents(pressEvInit));
            EEG_table.Pressure(2) = extractBetweenEvents(EEG, uniqueEvents(pressEvFin));
            EEG_table.Pressure(3) = extractBetweenEvents(EEG, uniqueEvents(pressCpmEv));

            % -------------------------------------------------

            cd(svloc);
            temp_file = [id,' -- preprocessed.mat'];
            save(temp_file, 'EEG', 'EEG_table', 'rejected', 'pcaComp', 'pcaScr', 'pcaPexp', 'rejectedPCA'); 
            
        end
end

%% helper functions 
function evs = uniqueEvents(evs)
    idx1 = 1; idx2 = idx1 + 1;
    while idx1 < size(evs,1)
        ev1 = evs(idx1,:); ev2 = evs(idx2,:);
        if (ev2(1).init_time <= ev1(2).init_time) & (ev2(2).init_time >= ev1(1).init_time)
            % events are not unique (ev2 is contained in or equal to ev1) 
            evs = evs([1:(idx2-1),(idx2+1):end],:);
        elseif idx2 < size(evs,1)
            idx2 = idx2 + 1;
        else
            idx1 = idx1 + 1; idx2 = idx1 + 1;
        end
    end
end

function outEEG = extractBetweenEvents(inEEG, evs)
    buf = .5; % s
    outEEG = repmat(inEEG, 1, size(evs,1));
    for idx = 1:size(evs,1)
        bound = [evs(idx,:).init_time];
        bound(1) = max(inEEG.xmin, bound(1)-buf);
        bound(2) = min(inEEG.xmax, bound(2)+buf);
        outEEG(idx) = pop_select(inEEG, 'time', bound);
    end
    outEEG = {outEEG};
end

function [startEvs, endEvs] = getTrainOfEvs(evs, t_inter_ev, t_before_ev)
    if nargin < 3
        t_before_ev = [];
    end
    isInitialBeforeTrain = length(t_before_ev);

    init_time = [evs.init_time];
    timeDiff = diff(init_time);
    intvlFromPrev = [inf, timeDiff]; intvlToNext = [timeDiff, inf];

    firstOfTrain = (intvlToNext < t_inter_ev) & (intvlFromPrev >= t_inter_ev);
    lastOfTrain = (intvlFromPrev < t_inter_ev) & (intvlToNext >= t_inter_ev);
    if isInitialBeforeTrain
        evBefore = (intvlToNext >= t_inter_ev) & (intvlToNext < t_before_ev) & (intvlFromPrev >= t_inter_ev);
    end

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
    if isInitialBeforeTrain
        for idx = 1:length(evStart)
            e0 = find(evBefore);
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
    end

    if isInitialBeforeTrain
        startEvs = ev0;
        startEvs(ev0 < 0) = evStart(ev0 < 0);
    else
        startEvs = evStart;
    end
    startEvs = evs(startEvs); endEvs = evs(evEnd);
end

function [evs, evsInit, evsFin] = getBaselineEvs(startEv, endEv, excludeTypes, splitTypes)
    evs = matchStartEndEvs(startEv, endEv);

    % handle events between start/end 
    toExclude = false(size(evs,1),1);
    for idx = 1:size(evs,1)
        evBetween = EEG.event( ...
            (eventTime > evs(idx,1).init_time ) & ...
            (eventTime < evs(idx,2).init_time) );
        toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
    end
    evs = evs(~toExclude,:);

    [evsInit, evsFin] = splitBeforeAfter(evs, splitTypes);
end

function evs = matchStartEndEvs(startEv, endEv)
    % order by time 
    if length(startEv) > 1
        [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
    end
    if length(endEv) > 1
        [~,ord] = sort([endEv.init_time]); endEv = endEv(ord);
    end
    % match start/end 
    evs = [];
    for startEv_ = startEv
        endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
        if ~isempty(endEvOpts)
            [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
            endEv_ = endEvOpts(endEv_);
            evs = [evs; startEv_, endEv_];
        end
    end
end

function [evsInit, evsFin] = splitBeforeAfter(evs, splitTypes)
    % designate before/after experiment 
    toExclude = false(size(evs,1),1);
    for idx = 1:size(evs,1)
        evBetween = EEG.event( ...
            (eventTime < evs(idx,1).init_time ) );
        toExclude(idx) = sum( arrayfun(@(e) sum(e.type == splitTypes), evBetween) );
    end
    evsInit = evs(~toExclude,:); evsFin = evs(toExclude,:);
end