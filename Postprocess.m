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
%eeglab

%% 
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
                [~,ord] = sort(dur);
                EEG_table{r,c} = {cur(ord)};
            end
        end
    end
    %}

    % CHANGES TO MAKE: make table with fourier transform struct for ALL
    % things in the original table 
    % then: PAF 

    % Inspect Baseline Spectra 
    selVars = {'BaselineOpen', 'BaselineClosed', 'BaselineIce'};
    selRows = {'before experiment', 'after experiment'};
    before_after_spectra(EEG_table, selVars, selRows, fn);
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

function [Ptbl, fig1, ax, fig2] = before_after_spectra(tbl, vars, rows, sttl, comparRows)
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

    fig1 = figure; 
    sgtitle(sttl);
    subtbl = tbl(ismember(tbl.Properties.RowNames, rows), ...
                 ismember(tbl.Properties.VariableNames, vars));
    nchan = 0; chlocs = [];

    vartypes = arrayfun(@(i) 'cell', 1:size(subtbl,2), 'UniformOutput', false);
    Ptbl = table('Size',size(subtbl),'VariableTypes',vartypes);
    Ptbl.Properties.VariableNames = subtbl.Properties.VariableNames;
    Ptbl.Properties.RowNames = subtbl.Properties.RowNames;

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
                [~,~,P,w] = fftPlot(eeg.data, eeg.srate);
                Ptbl{r,v} = {P, w};
    
                semilogy(w, P); grid on;
                xlabel('\omega (Hz)'); ylabel('log[Pwr (\muV^2 s^2)]');
            end
            title([vname,' ',rttl]);

            idx = idx + 1;
        end
    end
    linkaxes(ax, 'xy');

    fig2 = figure;
    idx = 1;
    for v = 1:H
        w_v = sort(unique([Ptbl{:,v}{:,2}]));
        P_v = zeros(nchan,length(w_v), W);
        for r = 1:W
            subplot(H,W+1,idx);
            tblItem = Ptbl(r,v); 
            rttl = tblItem.Properties.RowNames{1};
            vname = tblItem.Properties.VariableNames{1};

            Pw = tblItem{1,1};
            if ~isempty(Pw)
                P = Pw{1}; w = Pw{2};
                P_v(:,:,r) = interp1(w, P, w_v);
                rho = corr(P');
                heatmap(rho);
            end
            title([vname,' ',rttl]);
            
            idx = idx + 1;
        end
        subplot(H,W+1,idx);
        title([rows{comparRows(1)},' vs ',rows{comparRows(2)}]);
        P_v = P_v(:,:,comparRows);
        rho = diag(corr(P_v(:,:,1)',P_v(:,:,2)'));
        topoplot(rho, chlocs);
        idx = idx + 1;
    end
end