%% Start eeglab
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Set Filepaths
clear 

home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)
load("BrainwaveFrequencyTable.mat");
global BandTableHz

[fn, fp] = uigetfile('*postprocessed.mat'); 
load([fp, '/', fn]);

%% first plot
meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);
plotOpts = {'Peak Freq (Hz)', 'Band Mean (\muV^2 s^2)', 'Band Density'};
fcnOpts = {@(w,P,band) peakFreq(w,P,band), meanAmp, ampDensity};
plotSel = listdlg_selectWrapper(plotOpts, 'single', 'Plot What?');
if sum(plotSel == 1:3)
    [bnd,bndname] = pickFrequency();
    yname = [bndname,' ',plotOpts{plotSel}];
    fcn0 = fcnOpts{plotSel};
    if sum(plotSel == [2,3])
        fcn = @(Spec, EEG) frqFcnEpoch(Spec, EEG, @(w,P) fcn0(w,P,bnd));
        if plotsel == 3
            % Density 
            ylims = [0 1];
        else
            ylims = [];
        end
    elseif plotSel == 1
        % peak freq
        fcn = @(Spec, EEG) peakFreqEpoch(Spec, EEG, bnd);
        if isa(bnd,'char') | isa(bnd,'string')
            ylims = band2freqs(bnd, BandTableHz);
        else
            ylims = bnd;
        end
    end
end
AllPlot_table = plotTbl(Epoch_table, EpochSpec_table, fcn, ylims, ...
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
%%
[TestPlot_table, Plot_table] = ...
    testTbl(AllPlot_table, BL, Epoch_table, {fn, yname}, testVars, testRows);

%% helper functions 

function subtbl = makeSubtbl(tbl, vars, rows)
    subtbl = tbl(ismember(tbl.Properties.RowNames, rows), ...
                 ismember(tbl.Properties.VariableNames, vars));
end

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

%% key functions 

function [tblOut, toTestTbl, epochsTbl, fig] = ...
    testTbl(toTestTbl, baseline, epochsTbl, sttl, vars, rows)
    if nargin > 4
        toTestTbl = makeSubtbl(toTestTbl, vars, rows);
        epochsTbl = makeSubtbl(epochsTbl, vars, rows);
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
                curEEG = epochsTbl{r,c}{1}(1);
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
                plotWithEvents(curTime, curTestOut, curEEG, [], ttl, 'p');

            end
            idx = idx + 1;
        end
    end
end

function [tblOut, epochTbl, epochSpecTbl, fig] = ...
    plotTbl(epochTbl, epochSpecTbl, fcn, ybound, yname, sttl, vars, rows)
    if nargin > 6
        epochTbl     = makeSubtbl(epochTbl,     vars, rows);
        epochSpecTbl = makeSubtbl(epochSpecTbl, vars, rows);
    elseif nargin < 6
        sttl = '';
        if nargin < 5
            yname = '';
            if nargin < 4
                ybound = [];
            end
        end
    end
    
    tblOut = epochTbl; fig = figure; sgtitle(sttl);
    W = height(epochTbl); H = width(epochTbl); idx = 1;
    for c = 1:H
        for r = 1:W
            curSpec = epochSpecTbl{r,c}{:};
            curEpoc = epochTbl{r,c}{:};
            if length(curEpoc) > 1
                EEG = curEpoc(1);
                [t,Y] = fcn(curSpec, curEpoc);
                tblOut{r,c} = {cat(3,t,Y)};
                subplot(H,W,idx); 
                ttl = [epochTbl.Properties.VariableNames{c},' ',epochTbl.Properties.RowNames{r}];
                plotWithEvents(t, Y, EEG, ybound, ttl, yname);
            elseif ~isempty(curEpoc)
                tblOut{r,c} = {curEpoc(1),[],[]};
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

    plt = plot(t,Y); hold on;
    title(ttl); xlabel('time (s)'); ylabel(ylbl);
    for ev = event
        lbl_sw = false;
        if ~isempty(ev.latency) & ~strcmp(ev.type,'boundary')
            initTime = ev.latency/srate;
            plot(initTime+[0,0], ybound, 'r', 'LineWidth',1.5);
            text(initTime, ybound(lbl_sw+1), ev.type);
            lbl_sw = ~lbl_sw;
        end
    end
end

function [times,Y] = frqFcnEpoch(epoch_Spec, epoch_EEG, fcn)
    Y = cell2mat( arrayfun(@(eegSpec) ...
        fcn(eegSpec.frequency1side, eegSpec.powerSpectrum), ...
        epoch_Spec(2:end), 'UniformOutput', false) );
    times = arrayfun(@(eeg) mean([eeg.xmin, eeg.xmax]), epoch_EEG(2:end));
    times = times'; Y = Y';
end

function [times,Freq] = peakFreqEpoch(epoch_Spec, epoch_EEG, bnd)
% inspect peakFreq across epochs (time)
    % index 1 = original 
    nchan = epoch_Spec(1).nbchan;
    PAFs = zeros(nchan, length(epoch_Spec)-1, 2);
    for chan = 1:nchan
        PAFs(chan,:,1) = arrayfun(@(eegSpec) ...
            peakFreq(eegSpec.frequency1side, eegSpec.powerSpectrum(chan,:), bnd), ...
            epoch_Spec(2:end));
        PAFs(chan,:,2) = arrayfun(@(eeg) ...
            mean([eeg.xmin, eeg.xmax]), epoch_EEG(2:end));
    end
    times = PAFs(:,:,2)'; Freq = PAFs(:,:,1)';
end