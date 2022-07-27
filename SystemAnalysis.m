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

%%
fp = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG/Data_Chronic Pain/Preprocessed AllCopies 2022-07-03/Postprocessed 2022-07-03 21.05.29 -- 4s epoch, 3.75s overlap/';
H01 = load([fp,'2021-11-11 CP H01 --- 20211111_1143.mat -- preprocessed.mat -- postprocessed.mat']);
P04 = load([fp,'2021-10-14 CP P04 --- 20211014_1308.mat -- preprocessed.mat -- postprocessed.mat']);

%%

H01_Epoch_t_PP    = H01.Epoch_table.PinPrick('after experiment');
H01_Epoch_t_PPCPM = H01.Epoch_table.PinPrick('CPM');
P04_Epoch_t_PP    = P04.Epoch_table.PinPrick('before experiment');
P04_Epoch_t_PPCPM = P04.Epoch_table.PinPrick('CPM');

H01_Epoch_w_PP    = H01.EpochSpec_table.PinPrick('after experiment');
H01_Epoch_w_PPCPM = H01.EpochSpec_table.PinPrick('CPM');
P04_Epoch_w_PP    = P04.EpochSpec_table.PinPrick('before experiment');
P04_Epoch_w_PPCPM = P04.EpochSpec_table.PinPrick('CPM');

H01_t_PP    = H01.EEG_table.PinPrick('after experiment');
H01_t_PPCPM = H01.EEG_table.PinPrick('CPM');
P04_t_PP    = P04.EEG_table.PinPrick('before experiment');
P04_t_PPCPM = P04.EEG_table.PinPrick('CPM');

H01_w_PP    = H01.Spec_table.PinPrick('after experiment');
H01_w_PPCPM = H01.Spec_table.PinPrick('CPM');
P04_w_PP    = P04.Spec_table.PinPrick('before experiment');
P04_w_PPCPM = P04.Spec_table.PinPrick('CPM');

H01_Epoch_t_BL = H01.Epoch_table.BaselineOpen('before experiment');
P04_Epoch_t_BL = P04.Epoch_table.BaselineOpen('before experiment');
H01_Epoch_w_BL = H01.EpochSpec_table.BaselineOpen('before experiment');
P04_Epoch_w_BL = P04.EpochSpec_table.BaselineOpen('before experiment');
H01_t_BL       = H01.EEG_table.BaselineOpen('before experiment');
P04_t_BL       = P04.EEG_table.BaselineOpen('before experiment');
H01_w_BL       = H01.Spec_table.BaselineOpen('before experiment');
P04_w_BL       = P04.Spec_table.BaselineOpen('before experiment');


H01_Epoch_t_PP    = H01_Epoch_t_PP{1};
H01_Epoch_t_PPCPM = H01_Epoch_t_PPCPM{1};
P04_Epoch_t_PP    = P04_Epoch_t_PP{1};
P04_Epoch_t_PPCPM = P04_Epoch_t_PPCPM{1};
H01_Epoch_w_PP    = H01_Epoch_w_PP{1};
H01_Epoch_w_PPCPM = H01_Epoch_w_PPCPM{1};
P04_Epoch_w_PP    = P04_Epoch_w_PP{1};
P04_Epoch_w_PPCPM = P04_Epoch_w_PPCPM{1};
H01_t_PP    = H01_t_PP{1};
H01_t_PPCPM = H01_t_PPCPM{1};
P04_t_PP    = P04_t_PP{1};
P04_t_PPCPM = P04_t_PPCPM{1};
H01_w_PP    = H01_w_PP{1};
H01_w_PPCPM = H01_w_PPCPM{1};
P04_w_PP    = P04_w_PP{1};
P04_w_PPCPM = P04_w_PPCPM{1};
H01_Epoch_t_BL = H01_Epoch_t_BL{1};
P04_Epoch_t_BL = P04_Epoch_t_BL{1};
H01_Epoch_w_BL = H01_Epoch_w_BL{1};
P04_Epoch_w_BL = P04_Epoch_w_BL{1};
H01_t_BL = H01_t_BL{1};
P04_t_BL = P04_t_BL{1};
H01_w_BL = H01_w_BL{1};
P04_w_BL = P04_w_BL{1};


H01_event_PP = H01_t_PP(1).event; 
H01_event_PP = H01_event_PP(strcmp({H01_event_PP.type}, '11'));
P04_event_PP = P04_t_PP(1).event; 
P04_event_PP = P04_event_PP(strcmp({P04_event_PP.type}, '11'));

H01_event_PPCPM = H01_t_PPCPM(1).event; 
H01_event_PPCPM = H01_event_PPCPM(strcmp({H01_event_PPCPM.type}, '10'));
P04_event_PPCPM = P04_t_PPCPM(1).event; 
P04_event_PPCPM = P04_event_PPCPM(strcmp({P04_event_PPCPM.type}, '10'));

H01_boundTimes_PP    = eventBoundTimes(H01_event_PP);
H01_boundTimes_PPCPM = eventBoundTimes(H01_event_PPCPM);
P04_boundTimes_PP    = eventBoundTimes(P04_event_PP);
P04_boundTimes_PPCPM = eventBoundTimes(P04_event_PPCPM);


clear H01 P04

%% select what to plot
meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);
plotOpts = {'Peak Freq (Hz)', 'Band Mean (\muV^2 s^2)', 'Band Density', ...
            'Node Degree', 'Connectivity Strength', 'Frequency Difference', ...
            'Frequency Assortativity', ...
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
if sum(plotSel == [4,7,8])
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
if sum(plotSel == 1:7)
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
    elseif plotSel == 7
        % assortativity 
        fcn{idx} = @(Spec, EEG) assortativity(Spec, EEG, bnd, Pcut, BandTableHz);
        ylims = [-1, 1];
    elseif plotSel == 4
        % node degree
        fcn{idx} = @(Spec, EEG) nodeDegree(Spec, EEG, Pcut, bnd, BandTableHz);
        ylims = 'numchan';
    elseif plotSel == 5
        % conn strength
        fcn{idx} = @(Spec, EEG) avg_node(Spec, EEG, ...
            @(SO) WPLI(SO.frequencySpectrum, [], bnd, SO.frequency2side, BandTableHz) );
        ylims = [];
    elseif plotSel == 6
        % Freq Diff 
        fcn{idx} = @(Spec, EEG) avg_node(Spec, EEG, ...
            @(SO) diffFreq(SO.frequency1side, SO.powerSpectrum, bnd, BandTableHz) );
        ylims = [0, 1];
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

%%
[Y,t] = fcn(H01_Epoch_w_PP, H01_Epoch_t_PP); t = t(:,1);
tt = [H01_event_PP([31,32,41]).latency]/500;
hbnd = tt(1:2) + [0,-.5]; ybnd = tt(2:3);
h_idx = (t <= hbnd(2)) & (t >= hbnd(1));
y_idx = (t <= ybnd(2)) & (t >= ybnd(1));
Yh = Y(h_idx,:); Yy = Y(y_idx,:);
th = t(h_idx);   ty = t(y_idx);

ypred = zeros(size(Yy));
tt = [H01_event_PP(32:41).latency]/500;
tt = tt - ty(1); ty = ty - ty(1); th = th - th(1);
figure; hold on;
for tDelta = tt
    tShift = th - tDelta;
    hShift = cell2mat( arrayfun(@(c) ...
            interp1(th, Yh(:,c), tShift, 'nearest', 'extrap'), ...
            1:size(Yh,2), 'UniformOutput',false) );
    plot(th, hShift);
    hDelta = cell2mat( arrayfun(@(c) ...
            interp1(th, hShift(:,c), ty, 'nearest', 'extrap'), ...
            1:size(Yh,2), 'UniformOutput',false) );
    ypred = ypred + hDelta;
end

figure; plot(ypred, Yy, '.');

%% helper functions

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

function T = eventBoundTimes(evs)
    intvl = diff([evs.init_time]); 
    intvlFromPrev = [inf, intvl]; intvlToNext = [intvl, inf];
    %inTrain = (intvlFromPrev < 2) | (intvlToNext < 2);
    %firstOfTrain = (intvl >= 2) & (intvl < 10);
    %prickBefore = intvl >= 10;
    firstOfTrain = (intvlToNext < 2) & (intvlFromPrev >= 2);
    lastOfTrain = (intvlFromPrev < 2) & (intvlToNext >= 2);
    prickBefore = (intvlToNext >= 2) & (intvlToNext < 10) & (intvlFromPrev >= 2);

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

function [Y,t,Y1,Y2] = fcnCorr(fcn1, fcn2, var1, var2)
    [Y1,t1] = fcn1(var1, var2); [Y2,t2] = fcn2(var1, var2);
    Y1 = Y1'; Y2 = Y2'; t1 = t1'; t2 = t2';

    t = sort(unique([t1(:);t2(:)]));
    if length(t) > 1
        Y1 = cell2mat( arrayfun(@(c) ...
            interp1(t1(c,:), Y1(c,:), t, 'linear', 'extrap'), ...
            1:size(Y1,1), 'UniformOutput',false) )';
        Y2 = cell2mat( arrayfun(@(c) ...
            interp1(t2(c,:), Y2(c,:), t, 'linear', 'extrap'), ...
            1:size(Y2,1), 'UniformOutput',false) )';
    end

    Y = zeros(length(t),1);
    for s = 1:length(t)
        Y(s) = corr(Y1(:,s), Y2(:,s), 'Type', 'Spearman', 'Rows', 'complete');
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

function [Y,t] = rawData(EEGObj, chanSelect)
    if nargin < 2
        chanSelect = [];
    end
    if ~isempty(chanSelect)
        EEGObj = pop_select(EEGObj, 'channel', chanSelect);
    end
    Y = EEGObj.data; t = EEGObj.times/1000; 
    t = repmat(t, size(Y,1), 1); 
    t = t'; Y = Y';
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

function [A,t] = avg_node(SpectObj, EEGObj, nodeFcn)
    % each channel's average nodeFcn with all others 
    A = zeros(SpectObj(1).nbchan,length(SpectObj));
    for s = 1:length(SpectObj)
        W = nodeFcn(SpectObj(s));
        A(:,s) = mean(W, 'omitnan');
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(A,1), 1); 
    A = A'; t = t';
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