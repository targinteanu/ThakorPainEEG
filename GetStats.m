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
[~,scanfiles] = listdlg_selectWrapper({scanfiles.name},'multiple');

cd(home); addpath(postproDir);
svloc = [postproDir,'/Stats ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
mkdir(svloc); svloc = [svloc,'/'];

%% select what to plot

[fcn, yname, ylims] = MeasurementSelector();

[~,testVars] = listdlg_selectWrapper(...
    {'TempStim', 'PinPrick', 'Pressure', 'BaselineOpen', 'BaselineClosed', 'BaselineIce'}, ...
    'multiple', 'Specify Variable');

%% baseline and calculations  
%table2baseline = @(tbl) [tbl.BaselineOpen('before experiment'),...
%                         tbl.BaselineOpen('after experiment')];
table2baseline = @(tbl) tbl.BaselineOpen('before experiment');
p = .05; % uncertainty level 
timeBetweenEvents = 5; timeAfterLast = 3; % seconds 
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
    n = arrayfun(@(c) sum(~isnan(BLY(:,c))), 1:size(BLY,2));
    dof = n - 1; % deg of freedom 
    SE = std(BLY, 'omitnan') ./ sqrt(n); % standard error
    BL_CI = mean(BLY, 'omitnan') + [-1; 0; 1] .* tinv(p, dof).*SE;
    BLs{subj,1} = BL_CI; BLs{subj,2} = BL;
    clear BL_Epoch BL_EpochSpec t Y idx dof SE BLY BL BL_CI n

    % run calculations on desired variables 
    disp('calculating desired variables')
    [tY_table, EEG_table] = ...
        fcnTbl(EEG_table, Epoch_table, EpochSpec_table, fcn, testVars);

    dataTables{subj,1} = tY_table; dataTables{subj,2} = EEG_table;
    clear tY_table EEG_table Epoch_table EpochSpec_table Spec_table
end

%% graph baselines of each subj 
fig(2) = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(['Baseline ',yname]);
fig(1) = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(['Baseline ',yname]);
for subj = 1:size(BLs,1)
    fn = scanfiles{subj};
    %{
    strend = find(fn == '-'); strend = strend((diff(strend)==1)); strend = max(1, strend(1)-2);
    pname = fn(1:strend);
    pname = pname((end-2):end);
    %}
    pname = fn(1:3);

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
        
        figure(fig(2)); 
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
%maxTrialDur = maxTrialDur + timeAfterLast;

%% segment into trials 
trialTables = cell(size(dataTables));
for subj = 1:size(dataTables,1)
    fn = scanfiles{subj}

    EEG_table = dataTables{subj,2};
    %EEG_table = makeSubtbl(EEG_table, v);
    EEG_trial = EEG_table;
    tY_table  = dataTables{subj,1};
    %tY_table  = makeSubtbl(tY_table,  v);
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
                
                if strcmp(v, 'PinPrick')
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
    fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle([yname,' ',v{:}]);
    idx3 = 0; 
    for subj = 1:length(scanfiles)
        idx3incr = false;
        fn = scanfiles{subj}
        %{
        strend = find(fn == '-'); strend = strend((diff(strend)==1)); strend = max(1, strend(1)-2);
        ylbl = fn(1:strend);
        ylbl = ylbl((end-2):end);
        %}
        ylbl = fn(1:3);

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
    fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle([yname,' ',v{:}]);
    idx3 = 0;
    for subj = 1:length(scanfiles)
        idx3incr = false;
        fn = scanfiles{subj}
        %{
        strend = find(fn == '-'); strend = strend((diff(strend)==1)); strend = max(1, strend(1)-2);
        ylbl = fn(1:strend);
        ylbl = ylbl((end-2):end);
        %}
        ylbl = fn(1:3);

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
                            if (length(y) > 1) & (length(bl) > 1)
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
    fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(sttl); 
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