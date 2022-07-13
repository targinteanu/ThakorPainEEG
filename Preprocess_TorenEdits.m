%% Set Filepaths
clear

% functions for manipuliating dirs:
% get only subfolders
subfoldersof = @(d) d([d.isdir] & ...
    ~strcmp({d.name}, '.') & ...
    ~strcmp({d.name}, '..'));

home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG/Data_Chronic Pain';
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

addpath '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG/';
svloc = [home,'/Preprocessed ',datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];

%% Start eeglab
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Preprocessing in EEGLAB ---------------------------------------------
mkdir(svloc);

extractBetweenEvents = @(inEEG, evs) ...
    {arrayfun(@(idx) pop_select(inEEG, 'time', [evs(idx,:).init_time]), 1:size(evs,1))};

for subj = 1:size(datafolders,1)
    cd(home)
    disp(datafolders(subj,1).name)
    cd(datafolders(subj,1).name);
    clearvars infolder
    infolder = dir; infolder = subfoldersof(infolder);
    for d1 = 1:size(infolder,1)        
        disp(infolder(d1,1).name)
        cd(infolder(d1,1).name);

        clearvars datasets eventsets
        datasets = dir('*.cnt');
        eventsets = dir('*.ev2');

        % reorder dirs by dates 
        if size(datasets,1) > 1
            [~,ord] = sort([datasets.datenum]); datasets = datasets(ord);
        end
        if size(eventsets,1) > 1
            [~,ord] = sort([eventsets.datenum]); eventsets = eventsets(ord);
        end

        % handle mismatching number of data/event files 
        if size(datasets,1) ~= size(eventsets,1)
            warning('mismatching .cnt and .ev2 files')
        end
        N = min(size(datasets,1), size(eventsets,1));

        for d2 = 1:N
            clearvars EEG rejected EEG_table ...
                pcaComp pcaScr pcaPexp rejectedPCA ...
                endEv endEv_ startEv startEv_ endEvOpts eventType eventTime ord toExclude ...
                iceEv iceEvInit iceEvFin openEv openEvInit openEvFin closedEv closedEvInit closedEvFin ...
                pressEv pressEvInit pressEvFin pressCpmEv prickEv prickEvInit prickEvFin prickCpmEv ...
                tempEv tempEvInit tempEvFin 

            id = [datafolders(subj,1).name,' --- ',infolder(d1,1).name,' --- ',datasets(d2,1).name];

            % import EEG and EV2
            EEG = pop_loadcnt(datasets(d2,1).name, 'dataformat', 'auto', 'memmapfile', '');
            EEG = eeg_checkset( EEG );
            EEG = pop_chanedit(EEG, 'lookup', [eeglabpath,'/plugins/dipfit/standard_BEM/elec/standard_1005.elc']);
            EEG = eeg_checkset( EEG );
            EEG = pop_importev2(EEG, eventsets(d2,1).name);

            %%{

            % initial modification 
            if EEG.srate > 512
                EEG = pop_resample( EEG, 512);
            end
            EEG = pop_select( EEG,'nochannel',{'CB1' 'CB2' 'HEO' 'VEO' 'HEOG' 'VEOG'});
            EEG = pop_reref( EEG, []);
            EEG = pop_eegfiltnew(EEG, 2, 100, 826, 0, [], 1);

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


            % extract desired events --------------------------
            eventType = [EEG.event.type]; eventTime = [EEG.event.init_time];

            % baseline eyes open: 0 to 1
                startEv = EEG.event(eventType == 0); endEv = EEG.event(eventType == 1);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                if length(endEv) > 1
                    [~,ord] = sort([endEv.init_time]); endEv = endEv(ord);
                end
                % match start/end 
                openEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        openEv = [openEv; startEv_, endEv_];
                    end
                end
                % handle events between start/end 
                excludeTypes = 1:13;
                toExclude = false(size(openEv,1),1);
                for idx = 1:size(openEv,1)
                    evBetween = EEG.event( ...
                        (eventTime > openEv(idx,1).init_time ) & ...
                        (eventTime < openEv(idx,2).init_time) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                openEv = openEv(~toExclude,:);
                % designate before/after experiment 
                excludeTypes = 3:13;
                toExclude = false(size(openEv,1),1);
                for idx = 1:size(openEv,1)
                    evBetween = EEG.event( ...
                        (eventTime < openEv(idx,1).init_time ) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                openEvInit = openEv(~toExclude,:); openEvFin = openEv(toExclude,:);

            % baseline eyes closed: 1 to 2
                startEv = EEG.event(eventType == 1); endEv = EEG.event(eventType == 2);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                if length(endEv) > 1
                    [~,ord] = sort([endEv.init_time]); endEv = endEv(ord);
                end
                % match start/end 
                closedEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        closedEv = [closedEv; startEv_, endEv_];
                    end
                end
                % handle events between start/end 
                excludeTypes = [0, 2:13];
                toExclude = false(size(closedEv,1),1);
                for idx = 1:size(closedEv,1)
                    evBetween = EEG.event( ...
                        (eventTime > closedEv(idx,1).init_time ) & ...
                        (eventTime < closedEv(idx,2).init_time) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                closedEv = closedEv(~toExclude,:);
                % designate before/after experiment 
                excludeTypes = 3:13;
                toExclude = false(size(closedEv,1),1);
                for idx = 1:size(closedEv,1)
                    evBetween = EEG.event( ...
                        (eventTime < closedEv(idx,1).init_time ) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                closedEvInit = closedEv(~toExclude,:); closedEvFin = closedEv(toExclude,:);

            % hand in ice: 6 to 7
                startEv = EEG.event(eventType == 6); endEv = EEG.event(eventType == 7);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                if length(endEv) > 1
                    [~,ord] = sort([endEv.init_time]); endEv = endEv(ord);
                end
                % match start/end 
                iceEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        iceEv = [iceEv; startEv_, endEv_];
                    end
                end
                % handle events between start/end 
                excludeTypes = [6, 7];
                toExclude = false(size(iceEv,1),1);
                for idx = 1:size(iceEv,1)
                    evBetween = EEG.event( ...
                        (eventTime > iceEv(idx,1).init_time ) & ...
                        (eventTime < iceEv(idx,2).init_time) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                iceEv = iceEv(~toExclude,:);

            % temp stim: 3
                startEv = EEG.event(eventType == 3); 
                excludeTypes = [3, 5, 9, 12, 13];
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                % match start/end 
                tempEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        tempEv = [tempEv; startEv_, endEv_];
                    end
                end
                % designate before/after ice experiment 
                excludeTypes = [6, 7, 8:10, 12, 13];
                toExclude = false(size(tempEv,1),1);
                for idx = 1:size(tempEv,1)
                    evBetween = EEG.event( ...
                        (eventTime < tempEv(idx,1).init_time ) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                tempEvInit = tempEv(~toExclude,:); tempEvFin = tempEv(toExclude,:);

            % pin prick: 11
                startEv = EEG.event(eventType == 11); 
                excludeTypes = [11, 5, 9, 12, 13];
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                % match start/end 
                prickEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        prickEv = [prickEv; startEv_, endEv_];
                    end
                end
                % designate before/after ice experiment 
                excludeTypes = [6, 7, 8:10, 12, 13];
                toExclude = false(size(prickEv,1),1);
                for idx = 1:size(prickEv,1)
                    evBetween = EEG.event( ...
                        (eventTime < prickEv(idx,1).init_time ) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                prickEvInit = prickEv(~toExclude,:); prickEvFin = prickEv(toExclude,:);

            % CPM pin prick: 10
                startEv = EEG.event(eventType == 10); 
                excludeTypes = [10, 5, 9, 12, 13];
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                % match start/end 
                prickCpmEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        prickCpmEv = [prickCpmEv; startEv_, endEv_];
                    end
                end
                % select only those within iceEv????

            % pressure: 4 (to 5)
                startEv = EEG.event(eventType == 4); 
                excludeTypes = [4, 5, 9, 12, 13];
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                % match start/end 
                pressEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        pressEv = [pressEv; startEv_, endEv_];
                    end
                end
                % designate before/after ice experiment 
                excludeTypes = [6, 7, 8:10, 12, 13];
                toExclude = false(size(pressEv,1),1);
                for idx = 1:size(pressEv,1)
                    evBetween = EEG.event( ...
                        (eventTime < pressEv(idx,1).init_time ) );
                    toExclude(idx) = sum( arrayfun(@(e) sum(e.type == excludeTypes), evBetween) );
                end
                pressEvInit = pressEv(~toExclude,:); pressEvFin = pressEv(toExclude,:);

            % CPM pressure: 8 (to 9)
                startEv = EEG.event(eventType == 8); 
                excludeTypes = [8, 5, 9, 12, 13];
                toExclude = arrayfun(@(e) sum(e.type == excludeTypes), EEG.event) ;
                endEv = EEG.event(~toExclude);
                % order by time 
                if length(startEv) > 1
                    [~,ord] = sort([startEv.init_time]); startEv = startEv(ord);
                end
                % match start/end 
                pressCpmEv = [];
                for startEv_ = startEv
                    endEvOpts = endEv([endEv.init_time] >= startEv_.init_time);
                    if ~isempty(endEvOpts)
                        [~, endEv_] = min([endEvOpts.init_time] - startEv_.init_time);
                        endEv_ = endEvOpts(endEv_);
                        pressCpmEv = [pressCpmEv; startEv_, endEv_];
                    end
                end
                % select only those within iceEv????

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
            save(temp_file, 'EEG', 'EEG_table', 'rejected', ...
                'pcaComp', 'pcaScr', 'pcaPexp', 'rejectedPCA'); 
            
        end
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
        elseif idx2 <= size(evs,1)
            idx2 = idx2 + 1;
        else
            idx1 = idx1 + 1; idx2 = idx1 + 1;
        end
    end
end