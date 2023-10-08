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

    epochT = 4; % s
    epoch_dt = 4; % s

cd(home); addpath(preproDir);
svloc = [preproDir,'/Postprocessed ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS'),...
    ' -- ',num2str(epochT),'s epoch, ',num2str(epochT-epoch_dt),'s overlap'];

%% Start eeglab
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% 
mkdir(svloc);

for subj = 1:size(scanfiles,1)
    fn = scanfiles(subj,1).name
    load(fn);

    % EEG_table = makeSubtbl(EEG_table, {'PinPrick', 'BaselineOpen'});

    % order by longest duration 
    %%{
    disp('Reordering by duration')
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
    disp('Obtaining frequency spectra')
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

    % segment time series into _s epochs 
    % 
    %%{
    disp('Segmenting epochs')
    Epoch_table = EEG_table;
    EpochSpec_table = Epoch_table;
    for r = 1:height(Epoch_table)
        for c = 1:width(Epoch_table)
            disp([EEG_table.Properties.VariableNames{c},' ',EEG_table.Properties.RowNames{r}])
            curEEGlist = EEG_table{r,c}{:};
            EpocList = cell(size(curEEGlist));
            SpecList = cell(size(EpocList));
            if ~isempty(curEEGlist)
                for lstIdx = 1:length(curEEGlist)
                    eeg = curEEGlist(lstIdx);
                    if ~isempty(eeg)
                        t = (eeg.xmin):epoch_dt:((eeg.xmax)-epochT);
                        curEpochs = repmat(eeg, size(t));
                        curSpecs = repmat(Spec_table{r,c}{:}(1), size(t));
                        for idx = 1:length(t)
                            disp(['Epoch ',num2str(idx),' of ',num2str(length(t)),...
                                ' (',num2str(100*idx/length(t),3),'%)'])
                            curEpoch = pop_select(eeg, 'time', t(idx)+[0,epochT]);
                            curEpoch.xmin = curEpoch.xmin + t(idx);
                            curEpoch.xmax = curEpoch.xmax + t(idx);
                            curEpoch.times = curEpoch.times + t(idx)*1000;
                            curEpochs(idx) = curEpoch;
                            curSpec = fftPlot(curEpoch.data, curEpoch.srate);
                            curSpec.chanlocs = curEpoch.chanlocs; curSpec.nbchan = curEpoch.nbchan;
                            curSpecs(idx) = curSpec;
                        end
                        EpocList{lstIdx} = curEpochs;
                        SpecList{lstIdx} = curSpecs;
                    end
                end
                Epoch_table{r,c} = EpocList;
                EpochSpec_table{r,c} = SpecList;
            end
        end
    end
    clear curEEGlist eeg t curEpoch curEpochs EpocList curSpec curSpecs SpecList
    %}

    save([svloc,'/',fn,' -- ','postprocessed.mat'], ...
        'EEG_table', 'Epoch_table', 'Spec_table', 'EpochSpec_table');

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

function subtbl = makeSubtbl(tbl, vars)
    subtbl = tbl(:, ismember(tbl.Properties.VariableNames, vars));
end