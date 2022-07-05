%% Start eeglab
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Set Filepaths
clear 

home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)

[fn, fp] = uigetfile('*postprocessed.mat'); 
load([fp, '/', fn]);

%% main plotting 
meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);

%% helper functions 
function subtbl = makeSubtbl(tbl, vars, rows)
    subtbl = tbl(ismember(tbl.Properties.RowNames, rows), ...
                 ismember(tbl.Properties.VariableNames, vars));
end

%% key functions 

function [tblOut, fig] = plotTbl(epochTbl, epochSpecTbl, fcn, ybound, yname, sttl, vars, rows)
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
                tblOut{r,c} = {EEG,t,Y};
                subplot(H,W,idx); 
                ttl = [epochTbl.Properties.VariableNames{c},' ',epochTbl.Properties.RowNames{r}];
                plotWithEvents(t, Y, EEG.event, ybound, ttl, yname);
            elseif ~isempty(curEpoc)
                tblOut{r,c} = {curEpoc(1),[],[]};
            end
            idx = idx + 1;
        end
    end
end


function plt = plotWithEvents(t, Y, event, ybound, ttl, ylbl)
    if nargin < 6
        ylbl = '';
        if nargin < 5
            ttl = '';
            if nargin < 4
                ybound = [];
                if nargin < 3
                    event = [];
                end
            end
        end
    end
    if isempty(ybound)
        ybound = [min(Y(:)), max(Y(:))];
    end

    plt = plot(t,Y); hold on;
    title(ttl); xlabel('time (s)'); ylabel(ylbl);
    for ev = event
        lbl_sw = false;
        if ~isempty(ev.latency) & ~strcmp(ev.type,'boundary')
            initTime = ev.latency/cur1.srate;
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