%% Start eeglab
clear 
clear global 

eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Set Filepaths 

home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)
load("BrainwaveFrequencyTable.mat");
global BandTableHz

[fn, fp] = uigetfile('*postprocessed.mat'); 
load([fp, '/', fn]);
% pull out one EEG to use for settings 
EEG0 = EEG_table.BaselineOpen('before experiment'); EEG0 = EEG0{1}(1);

%% first plot
clear fcn fcn0 yname ylims idx Pcut

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
        ylims = [0,EEG0.nbchan];
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

AllPlot_table = plotTbl(EEG_table, Epoch_table, EpochSpec_table, fcn, ylims, ...
    yname, fn);

%% select baseline and compare 
[~,blVars] = listdlg_selectWrapper(Epoch_table.Properties.VariableNames, ...
                                    'multiple', 'Specify Baseline(s)');
[~,blRows] = listdlg_selectWrapper(Epoch_table.Properties.RowNames, ...
                                    'multiple', 'Specify Baseline(s)');
BL_table = makeSubtbl(AllPlot_table, blVars, blRows);
BL = BL_table{:,:};
BL = cell2mat(reshape(BL,[],1));

[~,testVars] = listdlg_selectWrapper(Epoch_table.Properties.VariableNames, ...
                                    'multiple', 'Specify Variable(s)');
[~,testRows] = listdlg_selectWrapper(Epoch_table.Properties.RowNames, ...
                                    'multiple', 'Specify Variable(s)');

[TestPlot_table, Plot_table] = ...
    testTbl(AllPlot_table, BL, EEG_table, {fn, yname}, testVars, testRows);

%% helper functions 

function subtbl = makeSubtbl(tbl, vars, rows)
    subtbl = tbl(ismember(tbl.Properties.RowNames, rows), ...
                 ismember(tbl.Properties.VariableNames, vars));
end

function t = getTimes(epoch_EEG)
    if isempty(epoch_EEG)
        t = [];
    else
        t = arrayfun(@(eeg) mean([eeg.xmin, eeg.xmax]), epoch_EEG);
    end
end

function [bandrange, bandname] = pickFrequency()
    global BandTableHz
    freqOpts = [BandTableHz.Properties.RowNames; 'custom'];
    [bandrange,ok] = listdlg('ListString',freqOpts, 'SelectionMode','single', ...
        'PromptString', 'Specify Frequency Band');
    while ~ok
        bandrange = questdlg('Use unbounded band?');
        ok = strcmp(bandrange, 'Yes');
        if ~ok
            [bandrange,ok] = listdlg('ListString',freqOpts, 'SelectionMode','single', ...
                'PromptString', 'Specify Frequency Band');
        else
            bandrange = 0;
        end
    end
    if bandrange == 0
        bandrange = []; bandname = '';
    elseif bandrange == length(freqOpts)
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
    if ~ok
        if strcmp(SelectionMode,'multiple')
            sel = questdlg('select all?');
            ok = strcmp(sel, 'Yes');
            if ~ok
                sel = []; 
            else
                sel = 1:length(list);
            end
        else
            sel = [];
        end
    end
    listOut = list(sel);
end

function [Y,t] = fcnCorr(fcn1, fcn2, var1, var2)
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

function [tblOut, eegTbl, epochTbl, epochSpecTbl, fig] = ...
    plotTbl(eegTbl, epochTbl, epochSpecTbl, fcn, ybound, yname, sttl, vars, rows)
    if nargin > 7
        eegTbl       = makeSubtbl(eegTbl,       vars, rows);
        epochTbl     = makeSubtbl(epochTbl,     vars, rows);
        epochSpecTbl = makeSubtbl(epochSpecTbl, vars, rows);
    elseif nargin < 7
        sttl = '';
        if nargin < 6
            yname = '';
            if nargin < 5
                ybound = [];
            end
        end
    end
    
    tblOut = epochTbl; fig = figure; sgtitle({sttl, yname});
    W = height(epochTbl); H = width(epochTbl); idx = 1;
    for c = 1:H
        for r = 1:W
            curSpec = epochSpecTbl{r,c}{:};
            curEpoc = epochTbl{r,c}{:};
            curEEG  = eegTbl{r,c}{:};
            if length(curEpoc) > 1
                EEG = curEEG(1);
                [Y,t] = fcn(curSpec, curEpoc);
                tblOut{r,c} = {cat(3,t,Y)};
                subplot(H,W,idx); 
                ttl = [epochTbl.Properties.VariableNames{c},' ',epochTbl.Properties.RowNames{r}];
                plotWithEvents(t, Y, EEG, ybound, ttl, '');
            elseif ~isempty(curEpoc)
                tblOut{r,c} = {[]};
            end
            idx = idx + 1;
        end
    end
end


function plt = plotWithEvents(t, Y, EEG, ybound, ttl, ylbl)
    if nargin < 6
        ylbl = '';
        if nargin < 5
            ttl = '';
            if nargin < 4
                ybound = [];
            end
        end
    end
    if isempty(ybound)
        ybound = [min(Y(:)), max(Y(:))];
    end
    if (nargin < 3) | isempty(EEG)
        event = []; srate = 1;
    else
        event = EEG.event; srate = EEG.srate;
    end

    plt = plot(t,Y,'.'); hold on; ylim(ybound);
    title(ttl); xlabel('time (s)'); ylabel(ylbl);
    lbl_sw = false;
    for ev = event
        if ~isempty(ev.latency) & ~strcmp(ev.type,'boundary')
            initTime = ev.latency/srate;
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