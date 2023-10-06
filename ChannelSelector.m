function [chanSelIdx, chanSelName, chanObjAll, chanObjSel] = ChannelSelector(scanfiles, DATATABLES)

% get all channel names of all EEGs 
chanObjAll = [];
for s = 1:length(scanfiles)
    sf = scanfiles{s};
    dataTables = DATATABLES{s};
    for subj = 1:length(sf)
        dataTable = dataTables{subj,1};
        for c = 1:width(dataTable)
            EEG = dataTable{1,c}{1};
            chanObjAll = [chanObjAll, EEG.chanlocs];
        end
    end
end
clear sf dataTable dataTables EEG

% eliminate repeats (case insensitive) 
[~,idx] = unique(upper({chanObjAll.labels})); 
chanObjAll = chanObjAll(idx);

% only include chans common to all EEGs 
idx = true(size(chanObjAll));
for s = 1:length(scanfiles)
    sf = scanfiles{s};
    dataTables = DATATABLES{s};
    for subj = 1:length(sf)
        dataTable = dataTables{subj,1};
        for c = 1:width(dataTable)
            EEG = dataTable{1,c}{1};
            for ch = 1:length(chanObjAll)
                idx(ch) = idx(ch) & ...
                    sum( strcmpi(chanObjAll(ch).labels, {EEG.chanlocs.labels}) );
            end
        end
    end
end
chanObjAll = chanObjAll(idx);
clear sf dataTable dataTables EEG idx

% select desired chans 
figure; hold on;
for chan = chanObjAll
    plot3(chan.X, chan.Y, chan.Z, '.', ...
        'Color', chanColor(chan, chanObjAll));
    text(chan.X, chan.Y, chan.Z, chan.labels, ...
        'Color', chanColor(chan, chanObjAll));
end
[chanSelIdx, chanSelName] = listdlg_selectWrapper({chanObjAll.labels}, ...
    'multiple', 'Select Channels:');
chanObjSel = chanObjAll(chanSelIdx);

%% helper functions 
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
    rgb = .7 * rgb;
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

end