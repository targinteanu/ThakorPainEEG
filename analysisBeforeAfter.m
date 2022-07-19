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
plotOpts = {'Channel Head Map', 'Channel Correlation', 'Channel Spectra', ...
            'Network Matrix', '3D Network'};
plotSel = listdlg_selectWrapper(plotOpts, 'multiple');
for ps = plotSel
    if ps == 1
        % pick function 
        fcnOptNames = {'Peak Frequency in Band', 'Mean Band Power', 'Band Relative Density', ...
                       'Node Degree', 'Connectivity Strength'};
        fcnSel = listdlg_selectWrapper(fcnOptNames, 'single'); 
        fcnSelName = fcnOptNames{fcnSel};

        % pick frequency band range
        [bnd,bndname] = pickFrequency;
        if sum(fcnSel == 1:3)
            meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
            ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);
            fcnOpts = {@(w,P,band) peakFreq(w,P,band), meanAmp, ampDensity};
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
        else
            if fcnSel == 4
                % node deg
                Pcut = inputdlg('Cutoff Percentile (%)','Network Cutoff Selection');
                Pcut = str2num(Pcut{1});
                fcnSel = @(SpectObj) nodeDegree(SpectObj, Pcut, bnd, BandTableHz);
                maplims = 'numchan';
            elseif fcnSel == 5
                % Conn Strength
                fcnSel = @(SpectObj) mean(...
                    WPLI(SpectObj.frequencySpectrum, [], bnd, SpectObj.frequency2side, BandTableHz), ...
                    'omitnan');
                maplims = 'maxmin';
            end

            network_headmap(Spec_table, {fn,fcnSelName}, fcnSel, maplims, ...
                selVars, selRows);
        end

    elseif ps == 2
        % corr
        before_after_corr(Spec_table, fn, selVars, selRows);

    elseif ps == 3
        % spect
        before_after_spectra(Spec_table, fn, selVars, selRows);

    elseif sum(ps == [4,5])
        % Network 
        show3Dmap = ps == 5;
        Mopts = {'WPLI', 'Binary Adjacency', 'Frequency Difference'};
        Msel = listdlg_selectWrapper(Mopts, 'multiple', 'Display What?');
        [bnd,bndname] = pickFrequency;
        for MS = Msel
            if sum(MS == [1,2])
                showBinary = MS == 2;
                if showBinary
                    Pcut = inputdlg('Cutoff Percentile (%)','Network Cutoff Selection');
                    Pcut = str2num(Pcut{1});
                else
                    Pcut = [];
                end
                before_after_WPLI(Spec_table, Pcut, bnd, BandTableHz, showBinary, show3Dmap, ...
                    {fn, [bndname,' WPLI']}, ...
                    selVars, selRows);
            elseif MS == 3
                before_after_freqDiff(Spec_table, bnd, show3Dmap, ...
                    {fn, [bndname,' band frequency difference']}, ...
                    selVars, selRows);
            end
        end
    end
end

%% helper functions 

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

function fig1 = network_headmap(tbl, sttl, fcnToMap, maplims, vars, rows)
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
                chlocs = eeg.chanlocs; nchan = eeg.nbchan;
                if strcmp(maplims, 'numchan')
                    ml = [0, nchan];
                else
                    ml = maplims;
                end
    
                topoplot(fcnToMap(eeg), chlocs, 'maplimits', ml, 'electrodes','labels'); colorbar;

            end
            title([vname,' ',rttl]);

            idx = idx + 1;
        end
    end
end

function nd = nodeDegree(SpectObj, cutoffPercentile, bnd, tbl)
    if nargin < 4
        tbl = [];
        if nargin < 3
            bnd = [];
            if nargin < 2
                cutoffPercentile = [];
            end
        end
    end
    [~,A] = WPLI(SpectObj.frequencySpectrum, cutoffPercentile, bnd, SpectObj.frequency2side, tbl);
    nd = sum(A);
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

function plot3Dnetwork(Hmp, chlocs)
    badch = arrayfun(@(ch) isempty(ch.X) | isempty(ch.Y) | isempty(ch.Z), chlocs);
    goodch = ~badch;
    chlocs = chlocs(goodch);
    Hmp = Hmp(goodch, goodch);
    for c1 = 1:(size(Hmp,1)-1)
        for c2 = (c1+1):size(Hmp,2)
            %lw = Hmp(c1,c2)*4.9/(max(Hmp(:))-min(Hmp(:))) + .1;
            %colr = ones(1,3) - Hmp(c1,c2)/max(Hmp(:));
            colr = Hmp(c1,c2)/max(Hmp(:));
            ch12 = chlocs([c1, c2]);
            %plot3([ch12.X],[ch12.Y],[ch12.Z],'LineWidth',.2,'Color',colr);
            ch12 = [ch12, fliplr(ch12)];
            patch([ch12.X] + [0,0,1,1], [ch12.Y] + [0,0,1,1], [ch12.Z] + [0,0,1,1], ...
                'k', 'FaceAlpha', colr, 'EdgeColor', 'none');
        end
    end
    text([chlocs.X],[chlocs.Y],[chlocs.Z],{chlocs.labels},'Color','b');
end

function fig = before_after_WPLI(tbl, Pcut, bnd, bndTbl, binaryOnly, show3Dmap, sttl, vars, rows)
    
        if nargin < 9
            subtbl = tbl;
            if nargin < 7
                sttl = '';
                if nargin < 6
                    show3Dmap = false;
                    if nargin < 5
                        binaryOnly = false;
                        if nargin < 4
                            bndTbl = [];
                            if nargin < 3
                                bnd = [];
                                if nargin < 2
                                    Pcut = [];
                                end
                            end
                        end
                    end
                end
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
                    [PLI, PLIA] = WPLI(F, Pcut, bnd, w, bndTbl);
                    PLIA = double(PLIA);
                    if binaryOnly
                        Hmp = PLIA;
                    else
                        Hmp = PLI;
                    end
                    chlocs = Fw.chanlocs;

                    if ~show3Dmap
                        heatmap({chlocs.labels}, {chlocs.labels}, Hmp);
                    else
                        hold on;
                        plot3Dnetwork(Hmp, chlocs);
                    end
                end
                title([vname,' ',rttl]);
            end
            
            idx = idx + 1;
        end
    end
end

function fig = before_after_freqDiff(tbl, bnd, show3Dmap, sttl, vars, rows)
    
        if nargin < 6
            subtbl = tbl;
            if nargin < 4
                sttl = '';
                if nargin < 3
                    show3Dmap = false;
                    if nargin < 2
                        bnd = [];
                    end
                end
            end
        else
            subtbl = makeSubtbl(tbl, vars, rows);
        end
        global BandTableHz

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
                    F = Fw.powerSpectrum; w = Fw.frequency1side;
                    chlocs = Fw.chanlocs; 

                    Y = diffFreq(w, F, bnd, BandTableHz);
                    if ~show3Dmap
                        heatmap({chlocs.labels}, {chlocs.labels}, Y);
                    else
                        hold on;
                        plot3Dnetwork(Y, chlocs);
                    end
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