%% Set Filepaths
clear

% functions for manipuliating dirs:
% get only subfolders
subfoldersof = @(d) d([d.isdir] & ...
    ~strcmp({d.name}, '.') & ...
    ~strcmp({d.name}, '..'));

home = 'C:\Users\targi\Desktop\Thakor Chronic Pain Data\Data_Chronic Pain';
cd(home)
datafolders = dir;
datafolders = subfoldersof(datafolders);

% user selects which scans to use
[sel, ok] = listdlg('ListString',{datafolders.name}, 'SelectionMode','multiple');
while ~ok
    sel = questdlg('select all?');
    ok = strcmp(sel, 'Yes');
    if ~ok
        [sel, ok] = listdlg('ListString',{datafolders.name}, 'SelectionMode','multiple');
    else
        sel = 1:length(datafolders);
    end
end
datafolders = datafolders(sel);

addpath 'C:\Users\targi\Documents\MATLAB\ThakorPainEEG';
svloc = [home,'/Preprocessed ',datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];

%% Start eeglab
eeglabpath = 'C:\Program Files\MATLAB\R2022a\eeglab2023.0';
addpath(eeglabpath)
eeglab

%% Preprocessing in EEGLAB ---------------------------------------------
mkdir(svloc);

for subj = 1:size(datafolders,1)
    cd(home)
    disp(datafolders(subj,1).name)
    cd(datafolders(subj,1).name);
    clearvars infolder
    infolder = dir; %infolder = subfoldersof(infolder);
    for d1 = 1:size(infolder,1)        
        disp(infolder(d1,1).name)
        cd(infolder(d1,1).name);

        clearvars datasets eventsets
        datasets = dir('*.cnt');
        eventsets = dir('*.ev2');

        % reorder dirs by dates 
        if size(datasets,1) > 1
            [~,ord] = sort([datasets.datenum]); datasets = datasets(ord);
        end
        if size(eventsets,1) > 1
            [~,ord] = sort([eventsets.datenum]); eventsets = eventsets(ord);
        end

        % handle mismatching number of data/event files 
        if size(datasets,1) ~= size(eventsets,1)
            warning('mismatching .cnt and .ev2 files')
        end
        N = min(size(datasets,1), size(eventsets,1));

        for d2 = 1:N
            clearvars EEG rejected EEG_table ...
                pcaComp pcaScr pcaPexp rejectedPCA ...
                endEv endEv_ startEv startEv_ endEvOpts eventType eventTime ord toExclude ...
                iceEv iceEvInit iceEvFin openEv openEvInit openEvFin closedEv closedEvInit closedEvFin ...
                pressEv pressEvInit pressEvFin pressCpmEv prickEv prickEvInit prickEvFin prickCpmEv ...
                tempEv tempEvInit tempEvFin 

            id = [datafolders(subj,1).name,' --- ',infolder(d1,1).name,' --- ',datasets(d2,1).name];

            % import EEG and EV2
            EEG = pop_loadcnt(datasets(d2,1).name, 'dataformat', 'auto', 'memmapfile', '');
            EEG = eeg_checkset( EEG );
            EEG = pop_importev2(EEG, eventsets(d2,1).name);

            %%
            cd(svloc);
            temp_file = [id,' -- preprocessed.mat'];
            save(temp_file, 'EEG'); 
            
        end
    end
end