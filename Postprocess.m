%% Set Filepaths
clear 

home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)
preproDir = uigetdir;

cd(preproDir);
scanfiles = dir; 
[sel, ok] = listdlg('ListString',{scanfiles.name}, 'SelectionMode','multiple');
while ~ok
    sel = questdlg('select all?');
    ok = strcmp(sel, 'Yes');
    if ~ok
        [sel, ok] = listdlg('ListString',{scanfiles.name}, 'SelectionMode','multiple');
    else
        sel = 1:length(scanfiles);
    end
end
scanfiles = scanfiles(sel);

cd(home); addpath(preproDir);
svloc = [preproDir,'/Postprocessed ',datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];

%% Start eeglab
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% 
meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);

for subj = 1:size(scanfiles,1)
    fn = scanfiles(subj,1).name
    load(fn);

    % order by longest duration 
    %%{
    for r = 1:height(EEG_table)
        for c = 1:width(EEG_table)
            cur = EEG_table{r,c}{:};
            if ~isempty(cur)
                dur = [cur.xmax] - [cur.xmin];
                [~,ord] = sort(dur); ord = fliplr(ord);
                EEG_table{r,c} = {cur(ord)};
            end
        end
    end
    clear cur dur ord
    %}

    % table of frequency spectra 
    %%{
    Spec_table = EEG_table;
    for r = 1:height(Spec_table)
        for c = 1:width(Spec_table)
            cur = Spec_table{r, c}{:};
            if ~isempty(cur)
                curSpecs = arrayfun(@(eeg) fftPlot(eeg.data, eeg.srate), cur);
                for idx = 1:length(curSpecs)
                    curSpecs(idx).chanlocs = cur(idx).chanlocs;
                    curSpecs(idx).nbchan = cur(idx).nbchan;
                end
                Spec_table{r,c} = {curSpecs};
            end
        end
    end
    clear cur curSpecs
    %}

    % Inspect Baseline Spectra 
    %%{
    selVars = {'BaselineOpen', 'BaselineClosed', 'BaselineIce'};
    selRows = {'before experiment', 'after experiment'};
    before_after_spectra(Spec_table, selVars, selRows, fn);
    %}

    % Inspect PAF 
    %%{
    selVars = {'BaselineOpen', 'BaselineClosed', 'BaselineIce'};
    selRows = {'before experiment', 'after experiment'};
    PAFheadmap(Spec_table, selVars, selRows, fn);
    %}

    % segment time series into 5s epochs 
    % (table of vectors where first element is original EEG and subsequent
    % are EEG epochs)
    %%{
    epochT = 4; % s
    Epoch_table = EEG_table;
    EpochSpec_table = Epoch_table;
    for r = 1:height(Epoch_table)
        for c = 1:width(Epoch_table)
            cur = EEG_table{r,c}{:};
            if ~isempty(cur)
                eeg = cur(1); % first / longest duration only
                t = (eeg.xmin):epochT:(eeg.xmax); 
                curEpochs = repmat(eeg, size(t)); 
                curSpecs = repmat(Spec_table{r,c}{:}(1), size(t));
                for idx = 2:length(t)
                    curEpoch = pop_select(eeg, 'time', t(idx-[1,0]));
                        curEpoch.xmin = curEpoch.xmin + t(idx-1);
                        curEpoch.xmax = curEpoch.xmax + t(idx-1);
                        curEpoch.times = curEpoch.times + t(idx-1)*1000;
                    curEpochs(idx) = curEpoch;
                    curSpec = fftPlot(curEpoch.data, curEpoch.srate);
                    curSpec.chanlocs = curEpoch.chanlocs; curSpec.nbchan = curEpoch.nbchan;
                    curSpecs(idx) = curSpec;
                end
                Epoch_table{r,c} = {curEpochs};
                EpochSpec_table{r,c} = {curSpecs};
            end
        end
    end
    clear cur eeg t curEpoch curEpochs curSpec curSpecs
    %}

    % inspect PAF across epochs (time) 
    %%{
    figure; sgtitle(fn);
    W = height(Epoch_table); H = width(Epoch_table); idx = 1;
    for c = 1:H
        for r = 1:W
            cur = Epoch_table{r,c}{:};
            if length(cur) > 1
                % cur(1) = original EEG
                nchan = cur(1).nbchan; 
                PAFs = zeros(nchan, length(cur)-1, 2);
                for chan = 1:nchan
                    PAFs(chan,:,1) = arrayfun(@(eegSpec) ...
                        getPAF(eegSpec.powerSpectrum(chan,:), ...
                        eegSpec.frequency1side), ...
                        EpochSpec_table{r,c}{:}(2:end));
                    PAFs(chan,:,2) = arrayfun(@(eeg) ...
                        mean([eeg.xmin, eeg.xmax]), Epoch_table{r,c}{:}(2:end));
                end

                subplot(H,W,idx); hold on;
                plot(PAFs(:,:,2)',PAFs(:,:,1)');
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('PAF (Hz)');
                for ev = cur(1).event
                    if ~isempty(ev.latency)
                        initTime = ev.latency/cur(1).srate;
                        plot(initTime+[0,.001], [9,11], 'r', 'LineWidth',1.5);
                        text(initTime, 9, ev.type);
                    end
                end
            end
            idx = idx + 1;
        end
    end
    clear idx H W PAFs chan nchan initTime ev cur
    %}

    % inspect mean beta/theta across epochs (time) 
    %%{
    fig(3) = figure; sgtitle(fn); 
    fig(2) = figure; sgtitle(fn);
    fig(1) = figure; sgtitle(fn);
    W = height(Epoch_table); H = width(Epoch_table); idx = 1;
    for c = 1:H
        for r = 1:W
            cur = EpochSpec_table{r,c}{:};
            if length(cur) > 1
                thetas = cell2mat( arrayfun(@(eegSpec) ...
                    meanAmp(eegSpec.frequency1side, eegSpec.powerSpectrum, [4,8]), ...
                    cur(2:end), 'UniformOutput', false) );
                betas = cell2mat( arrayfun(@(eegSpec) ...
                    meanAmp(eegSpec.frequency1side, eegSpec.powerSpectrum, [13,30]), ...
                    cur(2:end), 'UniformOutput', false) );
                alphas = cell2mat( arrayfun(@(eegSpec) ...
                    meanAmp(eegSpec.frequency1side, eegSpec.powerSpectrum, [9,11]), ...
                    cur(2:end), 'UniformOutput', false) );
                t = arrayfun(@(eeg) mean([eeg.xmin, eeg.xmax]), Epoch_table{r,c}{:}(2:end));

                cur1 = Epoch_table{r,c}{:}(1);
                figure(fig(1));
                subplot(H,W,idx); hold on;
                plot(t',thetas'); ymin = min(thetas(:)); ymax = max(thetas(:));
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('mean \theta (\muV^2 s^2)');
                    for ev = cur1.event
                        if ~isempty(ev.latency)
                            initTime = ev.latency/cur1.srate;
                            plot(initTime+[0,.001], [ymin,ymax], 'r', 'LineWidth',1.5);
                            text(initTime, ymin, ev.type);
                        end
                    end
                figure(fig(2));
                subplot(H,W,idx); hold on;
                plot(t',betas'); ymin = min(betas(:)); ymax = max(betas(:));
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('mean \beta (\muV^2 s^2)');
                    for ev = cur1.event
                        if ~isempty(ev.latency)
                            initTime = ev.latency/cur1.srate;
                            plot(initTime+[0,.001], [ymin,ymax], 'r', 'LineWidth',1.5);
                            text(initTime, ymin, ev.type);
                        end
                    end
                figure(fig(3));
                subplot(H,W,idx); hold on;
                %{
                Y = (thetas./betas); ymin = min(Y(:)); ymax = max(Y(:));
                plot(t',Y');
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('mean \theta / mean \beta');
                    for ev = cur1.event
                        if ~isempty(ev.latency)
                            initTime = ev.latency/cur1.srate;
                            plot(initTime+[0,.001], [ymin,ymax], 'r', 'LineWidth',1.5);
                            text(initTime, ymin, ev.type);
                        end
                    end
                %}
                plot(t',alphas'); ymin = min(alphas(:)); ymax = max(alphas(:));
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('mean \alpha (\muV^2 s^2)');
                    for ev = cur1.event
                        if ~isempty(ev.latency)
                            initTime = ev.latency/cur1.srate;
                            plot(initTime+[0,.001], [ymin,ymax], 'r', 'LineWidth',1.5);
                            text(initTime, ymin, ev.type);
                        end
                    end
            end
            idx = idx + 1;
        end
    end
    clear idx H W betas thetas alphas t ev f fig initTime cur cur1 Y ymin ymax
    %}

    fig(3) = figure; sgtitle(fn); 
    fig(2) = figure; sgtitle(fn);
    fig(1) = figure; sgtitle(fn);
    W = height(Epoch_table); H = width(Epoch_table); idx = 1;
    for c = 1:H
        for r = 1:W
            cur = EpochSpec_table{r,c}{:};
            if length(cur) > 1
                thetas = cell2mat( arrayfun(@(eegSpec) ...
                    ampDensity(eegSpec.frequency1side, eegSpec.powerSpectrum, [4,8]), ...
                    cur(2:end), 'UniformOutput', false) );
                betas = cell2mat( arrayfun(@(eegSpec) ...
                    ampDensity(eegSpec.frequency1side, eegSpec.powerSpectrum, [13,30]), ...
                    cur(2:end), 'UniformOutput', false) );
                alphas = cell2mat( arrayfun(@(eegSpec) ...
                    ampDensity(eegSpec.frequency1side, eegSpec.powerSpectrum, [9,11]), ...
                    cur(2:end), 'UniformOutput', false) );
                t = arrayfun(@(eeg) mean([eeg.xmin, eeg.xmax]), Epoch_table{r,c}{:}(2:end));

                cur1 = Epoch_table{r,c}{:}(1);
                figure(fig(1));
                subplot(H,W,idx); hold on;
                plot(t',thetas'); 
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('\theta density');
                    for ev = cur1.event
                        if ~isempty(ev.latency)
                            initTime = ev.latency/cur1.srate;
                            plot(initTime+[0,.001], [0,1], 'r', 'LineWidth',1.5);
                            text(initTime, 0, ev.type);
                        end
                    end
                figure(fig(2));
                subplot(H,W,idx); hold on;
                plot(t',betas'); 
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('\beta density');
                    for ev = cur1.event
                        if ~isempty(ev.latency)
                            initTime = ev.latency/cur1.srate;
                            plot(initTime+[0,.001], [0,1], 'r', 'LineWidth',1.5);
                            text(initTime, 0, ev.type);
                        end
                    end
                figure(fig(3));
                subplot(H,W,idx); hold on;
                plot(t',alphas'); 
                title([Epoch_table.Properties.VariableNames{c},' ',Epoch_table.Properties.RowNames{r}]);
                xlabel('time (s)'); ylabel('\alpha density');
                    for ev = cur1.event
                        if ~isempty(ev.latency)
                            initTime = ev.latency/cur1.srate;
                            plot(initTime+[0,.001], [0,1], 'r', 'LineWidth',1.5);
                            text(initTime, 0, ev.type);
                        end
                    end
            end
            idx = idx + 1;
        end
    end
    clear idx H W betas thetas alphas t ev f fig initTime cur cur1 Y

    clear EEG_table Epoch_table Spec_table EpochSpec_table
end

%% helper functions
function [SpecStruct, Y, w, P, wP] = fftPlot(y, fs)
    L = size(y,2);
    Y = fft(y');
    P = abs(Y).^2; P = P'; P = P(:, 1:L/2+1);
    wP = fs*(0:(L/2))/L;
    w = [fliplr(-wP(2:end)), wP];
    Y = fftshift(Y); Y = Y';
    SpecStruct.frequencySpectrum = Y; SpecStruct.frequency2side = w; 
    SpecStruct.powerSpectrum = P; SpecStruct.frequency1side = wP;
end

function [CoG, maxval, pklocs] = getPAF(P, w, alpha_bounds)
    if nargin < 3
        alpha_bounds = [9, 11]; % Hz
    end
    aband = (w >= alpha_bounds(1)) & (w <= alpha_bounds(2));
    P = P(aband); w = w(aband);
    [~,maxIdx] = max(P); maxval = w(maxIdx);
    [pkP,pklocs] = findpeaks(P, w);
    if length(pklocs) > 1
        % order by magnitude of peak
        [~,ord] = sort(pkP);
        pklocs = pklocs(fliplr(ord));
    end
    CoG = sum(w.*P)/sum(P);
end

function fig1 = PAFheadmap(tbl, vars, rows, sttl)
    fig1 = figure; sgtitle(sttl);
    subtbl = tbl(ismember(tbl.Properties.RowNames, rows), ...
                 ismember(tbl.Properties.VariableNames, vars));

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

                PAF = arrayfun(@(c) getPAF(P(c,:),w), 1:nchan);
    
                topoplot(PAF, chlocs, 'maplimits', [9, 11]); colorbar;
            end
            title([vname,' ',rttl]);

            idx = idx + 1;
        end
    end
end

function [fig1, ax, fig2] = before_after_spectra(tbl, vars, rows, sttl, comparRows)
    if nargin < 5
        comparRows = [1,2];
        if nargin < 4
            sttl = '';
            if nargin < 3
                rows = tbl.Properties.RowNames;
                if nargin < 2
                    vars = tbl.Properties.VariableNames;
                end
            end
        end
    end

    fig1 = figure; sgtitle(sttl);
    subtbl = tbl(ismember(tbl.Properties.RowNames, rows), ...
                 ismember(tbl.Properties.VariableNames, vars));
    nchan = 0; chlocs = [];

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
                nchan = eeg.nbchan; chlocs = eeg.chanlocs;
                P = eeg.powerSpectrum; w = eeg.frequency1side;
    
                semilogy(w, P); grid on;
                xlabel('\omega (Hz)'); ylabel('log[Pwr (\muV^2 s^2)]');
            end
            title([vname,' ',rttl]);

            idx = idx + 1;
        end
    end
    linkaxes(ax, 'xy');

    fig2 = figure; sgtitle(sttl);
    idx = 1;
    for v = 1:H
        tblRow = [subtbl{:,v}{:}]; 
        if ~isempty(tblRow)
            w_v = sort(unique([tblRow.frequency1side]));
        else
            w_v = [];
        end
        P_v = zeros(nchan,length(w_v), W);
        for r = 1:W
            subplot(H,W+1,idx);
            tblItem = subtbl(r,v); 
            rttl = tblItem.Properties.RowNames{1};
            vname = tblItem.Properties.VariableNames{1};

            Pw = tblItem{1,1}{:};
            if ~isempty(Pw)
                P = Pw.powerSpectrum; w = Pw.frequency1side;
                P_v(:,:,r) = interp1(w', P', w_v', 'linear', 'extrap')'; 
                rho = corr(P');
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
                topoplot(rho, chlocs, 'maplimits', [0, 1]); colorbar;
            end
        end
        idx = idx + 1;
    end
end