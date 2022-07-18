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

%% decide which things to show
varOpts = Spec_table.Properties.VariableNames; rowOpts = Spec_table.Properties.RowNames;
[~,selVars] = listdlg_selectWrapper(varOpts, 'multiple');
[~,selRows] = listdlg_selectWrapper(rowOpts, 'multiple');

%% main plotting 
plotOpts = {'Channel Head Map', 'Channel Correlation', 'Channel Spectra', 'WPLI'};
plotSel = listdlg_selectWrapper(plotOpts, 'multiple');
for ps = plotSel
    if ps == 1
        % pick frequency band range and function 
        [bnd,bndname] = pickFrequency;
        
        meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
        ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);
        fcnOpts = {@(w,P,band) peakFreq(w,P,band), meanAmp, ampDensity};
        fcnOptNames = {'Peak Frequency in Band', 'Mean Band Power', 'Band Relative Density'};
        fcnSel = listdlg_selectWrapper(fcnOptNames, 'single'); 
        fcnSelName = fcnOptNames{fcnSel};
        if fcnSel == 1
            % map limits based on band
            if isa(bnd,'char') | isa(bnd,'string')
                maplims = band2freqs(bnd, BandTableHz);
            else
                maplims = bnd;
            end
        elseif fcnSel == 3
            maplims = [0,1];
            if isa(bnd,'char') | isa(bnd,'string')
                bnd = band2freqs(bnd, BandTableHz);
            end
        elseif fcnSel == 2
            maplims = 'maxmin';
            if isa(bnd,'char') | isa(bnd,'string')
                bnd = band2freqs(bnd, BandTableHz);
            end
        end
        fcnSel = fcnOpts{fcnSel}; fcnSel = @(w,P) fcnSel(w,P,bnd);

        headmap(Spec_table, {fn,[bndname,' ',fcnSelName]}, ...
            fcnSel, maplims, selVars, selRows);

    elseif ps == 2
        % corr
        before_after_corr(Spec_table, fn, selVars, selRows);

    elseif ps == 3
        % spect
        before_after_spectra(Spec_table, fn, selVars, selRows);

    elseif ps == 4
        % WPLI
        before_after_WPLI(Spec_table, fn, selVars, selRows);

    end
end

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
        bandrange = freqOpts{bandrange}; bandname = bandrange;
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

function subtbl = makeSubtbl(tbl, vars, rows)
    subtbl = tbl(ismember(tbl.Properties.RowNames, rows), ...
                 ismember(tbl.Properties.VariableNames, vars));
end

%% key functions 

function fig1 = headmap(tbl, sttl, fcnToMap, maplims, vars, rows)
    fig1 = figure; sgtitle(sttl);

    if nargin > 4
        subtbl = makeSubtbl(tbl, vars, rows);
    else
        subtbl = tbl;
    end

    W = height(subtbl); H = width(subtbl);
    idx = 1;
    for v = 1:H
        for r = 1:W
            subplot(H,W,idx);

            tblItem = subtbl(r,v);
            rttl = tblItem.Properties.RowNames{1};
            vname = tblItem.Properties.VariableNames{1};
            eeg = tblItem{1,1}{1}; 
            if ~isempty(eeg)
                eeg = eeg(1); % first / longest duration only
                nchan = eeg.nbchan; chlocs = eeg.chanlocs;
                P = eeg.powerSpectrum; w = eeg.frequency1side;

                PAF = arrayfun(@(c) fcnToMap(w,P(c,:)), 1:nchan);
    
                topoplot(PAF, chlocs, 'maplimits', maplims, 'electrodes','labels'); colorbar;

            end
            title([vname,' ',rttl]);

            idx = idx + 1;
        end
    end
end

function fig2 = before_after_corr(tbl, sttl, vars, rows, comparRows)
    if nargin < 5
        comparRows = [1,2];
        if nargin < 4
            subtbl = tbl;
            if nargin < 2
                sttl = '';
            end
        else
            subtbl = makeSubtbl(tbl, vars, rows);
        end
    end

    fig2 = figure; sgtitle(sttl);
    idx = 1;
    W = height(subtbl); H = width(subtbl);
    for v = 1:H
        tblRow = [subtbl{:,v}{:}]; 
        if ~isempty(tblRow)
            w_v = sort(unique([tblRow.frequency1side]));
        else
            w_v = [];
        end
        P_v = zeros(min([tblRow.nbchan]),length(w_v), W);
        for r = 1:W
            subplot(H,W+1,idx);
            tblItem = subtbl(r,v); 
            rttl = tblItem.Properties.RowNames{1};
            vname = tblItem.Properties.VariableNames{1};

            chlocs = [];

            Pw = tblItem{1,1}{:};
            if ~isempty(Pw)
                P = Pw.powerSpectrum; w = Pw.frequency1side;
                P_v(:,:,r) = interp1(w', P', w_v', 'linear', 'extrap')'; 
                rho = corr(P');
                chlocs = Pw.chanlocs;
                heatmap({chlocs.labels}, {chlocs.labels}, rho);
            end
            title([vname,' ',rttl]);
            
            idx = idx + 1;
        end
        subplot(H,W+1,idx);
        title([rows{comparRows(1)},' vs ',rows{comparRows(2)}]);
        P_v = P_v(:,:,comparRows);
        if ~isempty(P_v)
            rho = diag(corr(P_v(:,:,1)',P_v(:,:,2)'))';
            %rho(isnan(rho)) = 0;
            if ~sum(isnan(rho))
                topoplot(rho, chlocs, 'maplimits', [-1, 1], 'electrodes','labels'); colorbar;
            end
        end
        idx = idx + 1;
    end
end

function fig = before_after_WPLI(tbl, sttl, vars, rows)
    
        if nargin < 4
            subtbl = tbl;
            if nargin < 2
                sttl = '';
            end
        else
            subtbl = makeSubtbl(tbl, vars, rows);
        end

    fig = figure; sgtitle(sttl);
    idx = 1;
    W = height(subtbl); H = width(subtbl);
    for v = 1:H
        for r = 1:W
            subplot(H,W,idx);
            tblItem = subtbl(r,v); 
            rttl = tblItem.Properties.RowNames{1};
            vname = tblItem.Properties.VariableNames{1};

            if ~isempty(tblItem{1,1})
                Fw = tblItem{1,1}{:};
                if ~isempty(Fw)
                    F = Fw.frequencySpectrum; w = Fw.frequency2side;
                    PLI = WPLI(F);
                    chlocs = Fw.chanlocs;
                    heatmap({chlocs.labels}, {chlocs.labels}, PLI);
                end
                title([vname,' ',rttl]);
            end
            
            idx = idx + 1;
        end
    end
end

function [fig1, ax] = before_after_spectra(tbl, sttl, vars, rows)
        if nargin < 4
            subtbl = tbl;
            if nargin < 2
                sttl = '';
            end
        else
            subtbl = makeSubtbl(tbl, vars, rows);
        end

    fig1 = figure; sgtitle(sttl);

    W = height(subtbl); H = width(subtbl);
    idx = 1;
    for v = 1:H
        for r = 1:W
            ax(idx) = subplot(H,W,idx);

            tblItem = subtbl(r,v);
            rttl = tblItem.Properties.RowNames{1};
            vname = tblItem.Properties.VariableNames{1};
            eeg = tblItem{1,1}{1}; 
            if ~isempty(eeg)
                eeg = eeg(1); % first / longest duration only
                P = eeg.powerSpectrum; w = eeg.frequency1side;
    
                semilogy(w, P); grid on;
                xlabel('\omega (Hz)'); ylabel('log[Pwr (\muV^2 s^2)]');
            end
            title([vname,' ',rttl]);

            idx = idx + 1;
        end
    end
    linkaxes(ax, 'xy');

end