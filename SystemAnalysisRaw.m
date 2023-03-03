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
svloc = [postproDir,'/System ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];

%% loading 
RPs = cell(length(scanfiles), 2, 2); Ys = cell(length(scanfiles), 2, 4);
ratings = cell(length(scanfiles), 2);
for subj = 1:length(scanfiles)
    fn = scanfiles{subj}
    PtLd = load(fn);

    tbls = {PtLd.Epoch_table; ...
            PtLd.EpochSpec_table; ...
            PtLd.EEG_table; ...
            PtLd.Spec_table};
    clear PtLd
    objStructs = cell(length(tbls), 3);
    % col 1 = no CPM; 2 = CPM; 3 = baseline

    %% without loading 

    % get PinPricks EEG/spectrum/epoch data -------------------------------
    for idx = 1:length(tbls)
        curTbl = tbls{idx};
        tempObj = curTbl.PinPrick('CPM');
        objStructs{idx,2} = tempObj{:}; 
        tempObj = curTbl.PinPrick('before experiment');
        tempObj1 = tempObj{:};
        if isempty(tempObj1)
            tempObj1 = [];
        end
        tempObj = curTbl.PinPrick('after experiment');
        tempObj = tempObj{:};
        if isempty(tempObj)
            tempObj = [];
        end
        objStructs{idx,1} = [tempObj1, tempObj];
        tempObj = curTbl.BaselineOpen('before experiment');
        objStructs{idx,3} = tempObj{:};
    end
    for cond = 1:3
        tempObj = objStructs{3,cond};
        if length(tempObj) > 1
            objStructs{3,cond} = pop_mergeset(tempObj(1), tempObj(2));
        end
    end
    clear tempObj tempObj1 

    % get events ----------------------------------------------------------
    evStructs = cell(3,2);
    tempEv = objStructs{3,1}; 
    if ~isempty(tempEv)
        tempEv = tempEv(1).event;
        evStructs{1,1} = tempEv(strcmp({tempEv.type}, '11'));
        evStructs{2,1} = eventBoundTimes(evStructs{1,1}); % <--- change to split evs by specific trial
        evStructs{3,1} = objStructs{3,1}(1).srate;
    end
    tempEv = objStructs{3,2}; 
    if ~isempty(tempEv)
        tempEv = tempEv(1).event;
        evStructs{1,2} = tempEv(strcmp({tempEv.type}, '10'));
        evStructs{2,2} = eventBoundTimes(evStructs{1,2});
        evStructs{3,2} = objStructs{3,2}(1).srate;
    end
    clear tempEv

    % calculations --------------------------------------------------------
    BL = objStructs{3,3}; BL_t = BL.times/1000'; BL = BL.data';
    BL_w = objStructs{4,3}; w = BL_w.frequency1side; BL_w = BL_w.powerSpectrum;

    for cond = 1:size(evStructs,2)
        % cond: 1 = no CPM, 2 = CPM
        if ~isempty(objStructs{1,cond})

        Y = objStructs{3,cond}; t = Y.times/1000'; Y = Y.data';
        Y_w = objStructs{4,cond}; w = Y_w.frequency1side; Y_w = Y_w.powerSpectrum;
        %Y = Y - mean(BL); % remove baseline
        %Y = (Y - mean(BL))./(std(BL)/sqrt(size(BL,1))); % t statistic
        chloc = objStructs{3,3}.chanlocs;

        srate = evStructs{3,cond};
        T = evStructs{2,cond};
        T = T(T(:,1)>=0, :);
        RP = zeros(2, size(Y,2), size(T,1)); RP(2,:,:) = 1;
        rtg = zeros(size(T,1), 3);
        Yin = nan(size(Y,1), size(Y,2), size(T,1)); Yout = Yin; Yir = Yin;
        X = zeros(size(Y,1), size(T,1));
        for trl = 1:size(T,1)
            trlEv = evStructs{1,cond}; % <-- change to get ev for specific trial 

            tsplit = T(trl,:);
            tt = [trlEv(tsplit).latency]/srate;

            hbnd = tt(1:2) + [0,-.15]; ybnd = tt(2:3);
            h_idx = (t <= hbnd(2)) & (t >= hbnd(1));
            y_idx = (t <= ybnd(2)) & (t >= ybnd(1));
            Yh = Y(h_idx,:); Yy = Y(y_idx,:);
            th = t(h_idx);   ty = t(y_idx);

            ypred = zeros(size(Yy));
            tt = [trlEv.latency]/srate;
            tt = tt( (tt <= ybnd(2)) & (tt >= ybnd(1)) );
            x = zeros(size(Yy,1),1);

            if ~isempty(ty) & (length(th)>1)
                tt = tt - ty(1); ty = ty - ty(1); th = th - th(1);

                for tDelta = tt
                    tShift = ty - tDelta;
                    hShift = cell2mat( arrayfun(@(c) ...
                        interp1(th, Yh(:,c), tShift, 'nearest', 'extrap')', ...
                        1:size(Yh,2), 'UniformOutput',false) );
                    ypred = ypred + hShift;

                    [~,xidx] = min(abs(tShift));
                    x(xidx) = 1;
                end

                for idx = 1:size(Yy,2)
                    [RP(1,idx,trl), RP(2,idx,trl)] = corr(ypred(:,idx), Yy(:,idx));
                end
                Yin(1:size(ypred,1),1:size(ypred,2),trl) = ypred; 
                Yout(1:size(Yy,1),1:size(Yy,2),trl) = Yy;
                Yir(1:size(Yh,1),1:size(Yh,2),trl) = Yh;
                X(1:length(x),trl) = x;
            end
        end
        RPs{subj, cond, 1} = RP; RPs{subj, cond, 2} = chloc;
        Ys{subj, cond, 1} = Yin; Ys{subj, cond, 2} = Yout; Ys{subj, cond, 3} = Yir; Ys{subj, cond, 4} = X;
        ratings{subj, cond} = rtg;
        clear RP chloc T Y t tt tsplit ypred Yy Yh th ty h_idx y_idx hbnd ybnd srate rtg 
        clear x X xidx hShift tShift tDelta Yir Yin Yout
        end
    end
    clear BL BLe

end

%% plotting 
maxNumTrl = 0;
for subj = 1:length(scanfiles)
    numTrl = size(RPs{subj, 1, 1},3) + size(RPs{subj, 2, 1},3);
    maxNumTrl = max(maxNumTrl, numTrl);
end
W = maxNumTrl; H = length(scanfiles);
fig(2) = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(['System p']);    % FIX
fig(1) = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(['System \rho']); % FIX
allchan = [];
for subj = 1:length(scanfiles)
    fn = scanfiles{subj};
    pname = fn(1:3);
    RP = RPs{subj, 1, 1};
    RP_CPM = RPs{subj, 2, 1};
    chloc = RPs{subj,1,2};
    chloc_CPM = RPs{subj,2,2};
    allchan = [allchan, chloc, chloc_CPM];

    if ~isempty(RP)
    for trl = 1:size(RP,3)
        idx = W*(subj-1) + trl;
        figure(fig(1)); subplot(H,W,idx);
        title([num2str(trl),' ',pname]);
        topoplot(RP(1,:,trl), chloc, 'maplimits',[-1 1]); colorbar;
        figure(fig(2)); subplot(H,W,idx);
        title([num2str(trl),' ',pname]);
        topoplot(RP(2,:,trl), chloc, 'maplimits',[0 1]); colorbar;
    end
    end

    if ~isempty(RP_CPM)
    for trl = 1:size(RP_CPM,3)
        idx = W*subj - trl + 1;
        figure(fig(1)); subplot(H,W,idx);
        title(['CPM ',num2str(trl),' ',pname]);
        topoplot(RP_CPM(1,:,trl), chloc, 'maplimits',[-1 1]); colorbar;
        figure(fig(2)); subplot(H,W,idx);
        title(['CPM ',num2str(trl),' ',pname]);
        topoplot(RP_CPM(2,:,trl), chloc, 'maplimits',[0 1]); colorbar;
    end
    end
end

%% channel selection 
figure; hold on;
[~,idx] = unique({allchan.labels}); 
allchan = allchan(idx);
for chan = allchan
    plot3(chan.X, chan.Y, chan.Z, '.', ...
        'Color', chanColor(chan, allchan));
    text(chan.X, chan.Y, chan.Z, chan.labels, ...
        'Color', chanColor(chan, allchan));
end
[chansel, chanselName] = listdlg_selectWrapper({allchan.labels}, ...
    'multiple', 'Select Channels:');

%% plotting in-out 
Mkr = {'o', 'x', '^', 's', 'v', 'p', '+', 'd', 'h', '*'};
fig(4) = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(['System IR']);     % FIX
fig(3) = figure('Units', 'Normalized', 'Position', [0 0 1 1]); sgtitle(['System in-out']); % FIX
H = 4; W = ceil(2*length(scanfiles)/H); % more robust?
for subj = 1:length(scanfiles)
    fn = scanfiles{subj};
    pname = fn(1:3);
    for cond = 1:2 % more robust?
        idx = mod(subj-1, W)+1 + (2*(ceil(subj/W)-1) + cond-1)*W;

        figure(fig(3)); ax3(idx) = subplot(H, W, idx); hold on;
        Yin = Ys{subj, cond, 1}; Yout = Ys{subj, cond, 2};
        chloc = RPs{subj, cond, 2};
        for trl = 1:size(Yin, 3)
            for ch = 1:size(Yin, 2)
                if sum(strcmp(chloc(ch).labels, chanselName))
                    plot(Yin(:,ch,trl), Yout(:,ch,trl), Mkr{trl}, ...
                        'Color', chanColor(chloc(ch), chloc));
                end
            end
        end
        grid on; 
        xlabel('LTI pred'); ylabel('actual');
        if cond == 1
            title(pname);
        elseif cond == 2
            title([pname,' CPM']);
        end

        figure(fig(4)); ax4(idx) = subplot(H,W,idx); hold on; 
        Yir = Ys{subj, cond, 3};
        for trl = 1:size(Yir, 3)
            for ch = 1:size(Yir, 2)
                if sum(strcmp(chloc(ch).labels, chanselName))
                    plot(Yir(:,ch,trl), ['-',Mkr{trl}], ...
                        'Color', chanColor(chloc(ch), chloc));
                end
            end
        end
        grid on; 
        xlabel('timestep'); % ylabel(yname); FIX
        if cond == 1
            title(pname);
        elseif cond == 2
            title([pname,' CPM']);
        end
    end
end
linkaxes(ax3); linkaxes(ax4);

%% saving 
mkdir(svloc); svloc = [svloc,'/'];
saveas(fig(1), [svloc,'correlation'], 'fig');
saveas(fig(2), [svloc,'p value'], 'fig');
saveas(fig(3), [svloc,'in-out'], 'fig');
saveas(fig(4), [svloc,'IR'], 'fig');

%% prep for sys iden 
subj = 4; cond = 1; ch = 3; trl = 1;
Yin = Ys{subj,cond,1}(:,ch,trl); Yout = Ys{subj,cond,2}(:,ch,trl);
Yh = Ys{subj,cond,3}(:,ch,trl);
X = Ys{subj,cond,4}(:,trl);
cens = isnan(Yin) | isnan(Yout);
Yin = Yin(~cens); Yout = Yout(~cens); X = X(~cens);

%% helper functions 
function T = eventBoundKNN(evs)
    % not functional - expand??

    intvl = diff([evs.init_time]);
    [idx, C] = kmeans(log(intvl)', 3);
    [~,ord] = sort(C); idx = ord(idx);

    T = 0;
end

function T = eventBoundTimes(evs)
    intvl = diff([evs.init_time]); 
    intvlFromPrev = [inf, intvl]; intvlToNext = [intvl, inf];
    
    tIS = 1.5; tIT = 10; % seconds

    firstOfTrain = (intvlToNext < tIS) & (intvlFromPrev >= tIS);
    lastOfTrain = (intvlFromPrev < tIS) & (intvlToNext >= tIS);
    prickBefore = (intvlToNext >= tIS) & (intvlToNext < tIT) & (intvlFromPrev >= 2);

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

    T = [ev0; evStart; evEnd]; T = T';
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