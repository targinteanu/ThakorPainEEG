%% Start eeglab
clear
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Set Filepaths
home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)

%%
fp = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG/Data_Chronic Pain/Preprocessed AllCopies 2022-07-03/Postprocessed 2022-07-03 21.05.29 -- 4s epoch, 3.75s overlap/';
H01 = load([fp,'2021-11-11 CP H01 --- 20211111_1143.mat -- preprocessed.mat -- postprocessed.mat']);
P04 = load([fp,'2021-10-14 CP P04 --- 20211014_1308.mat -- preprocessed.mat -- postprocessed.mat']);

%%

H01_Epoch_t_PP    = H01.Epoch_table.PinPrick('after experiment');
H01_Epoch_t_PPCPM = H01.Epoch_table.PinPrick('CPM');
P04_Epoch_t_PP    = P04.Epoch_table.PinPrick('before experiment');
P04_Epoch_t_PPCPM = P04.Epoch_table.PinPrick('CPM');

H01_Epoch_w_PP    = H01.EpochSpec_table.PinPrick('after experiment');
H01_Epoch_w_PPCPM = H01.EpochSpec_table.PinPrick('CPM');
P04_Epoch_w_PP    = P04.EpochSpec_table.PinPrick('before experiment');
P04_Epoch_w_PPCPM = P04.EpochSpec_table.PinPrick('CPM');

H01_t_PP    = H01.EEG_table.PinPrick('after experiment');
H01_t_PPCPM = H01.EEG_table.PinPrick('CPM');
P04_t_PP    = P04.EEG_table.PinPrick('before experiment');
P04_t_PPCPM = P04.EEG_table.PinPrick('CPM');

H01_w_PP    = H01.Spec_table.PinPrick('after experiment');
H01_w_PPCPM = H01.Spec_table.PinPrick('CPM');
P04_w_PP    = P04.Spec_table.PinPrick('before experiment');
P04_w_PPCPM = P04.Spec_table.PinPrick('CPM');

H01_Epoch_t_BL = H01.Epoch_table.BaselineOpen('before experiment');
P04_Epoch_t_BL = P04.Epoch_table.BaselineOpen('before experiment');
H01_Epoch_w_BL = H01.EpochSpec_table.BaselineOpen('before experiment');
P04_Epoch_w_BL = P04.EpochSpec_table.BaselineOpen('before experiment');
H01_t_BL       = H01.EEG_table.BaselineOpen('before experiment');
P04_t_BL       = P04.EEG_table.BaselineOpen('before experiment');
H01_w_BL       = H01.Spec_table.BaselineOpen('before experiment');
P04_w_BL       = P04.Spec_table.BaselineOpen('before experiment');


H01_Epoch_t_PP    = H01_Epoch_t_PP{1};
H01_Epoch_t_PPCPM = H01_Epoch_t_PPCPM{1};
P04_Epoch_t_PP    = P04_Epoch_t_PP{1};
P04_Epoch_t_PPCPM = P04_Epoch_t_PPCPM{1};
H01_Epoch_w_PP    = H01_Epoch_w_PP{1};
H01_Epoch_w_PPCPM = H01_Epoch_w_PPCPM{1};
P04_Epoch_w_PP    = P04_Epoch_w_PP{1};
P04_Epoch_w_PPCPM = P04_Epoch_w_PPCPM{1};
H01_t_PP    = H01_t_PP{1};
H01_t_PPCPM = H01_t_PPCPM{1};
P04_t_PP    = P04_t_PP{1};
P04_t_PPCPM = P04_t_PPCPM{1};
H01_w_PP    = H01_w_PP{1};
H01_w_PPCPM = H01_w_PPCPM{1};
P04_w_PP    = P04_w_PP{1};
P04_w_PPCPM = P04_w_PPCPM{1};
H01_Epoch_t_BL = H01_Epoch_t_BL{1};
P04_Epoch_t_BL = P04_Epoch_t_BL{1};
H01_Epoch_w_BL = H01_Epoch_w_BL{1};
P04_Epoch_w_BL = P04_Epoch_w_BL{1};
H01_t_BL = H01_t_BL{1};
P04_t_BL = P04_t_BL{1};
H01_w_BL = H01_w_BL{1};
P04_w_BL = P04_w_BL{1};


H01_event_PP = H01_t_PP(1).event; 
H01_event_PP = H01_event_PP(strcmp({H01_event_PP.type}, '11'));
P04_event_PP = P04_t_PP(1).event; 
P04_event_PP = P04_event_PP(strcmp({P04_event_PP.type}, '11'));

H01_event_PPCPM = H01_t_PPCPM(1).event; 
H01_event_PPCPM = H01_event_PPCPM(strcmp({H01_event_PPCPM.type}, '10'));
P04_event_PPCPM = P04_t_PPCPM(1).event; 
P04_event_PPCPM = P04_event_PPCPM(strcmp({P04_event_PPCPM.type}, '10'));

H01_boundTimes_PP    = eventBoundTimes(H01_event_PP);
H01_boundTimes_PPCPM = eventBoundTimes(H01_event_PPCPM);
P04_boundTimes_PP    = eventBoundTimes(P04_event_PP);
P04_boundTimes_PPCPM = eventBoundTimes(P04_event_PPCPM);


clear H01 P04

%% select what to plot
[fcn, yname, ylims] = MeasurementSelector();

%%
[Y,t] = fcn(H01_Epoch_w_PP, H01_Epoch_t_PP); t = t(:,1);
tt = [H01_event_PP([31,32,41]).latency]/500;
hbnd = tt(1:2) + [0,-.5]; ybnd = tt(2:3);
h_idx = (t <= hbnd(2)) & (t >= hbnd(1));
y_idx = (t <= ybnd(2)) & (t >= ybnd(1));
Yh = Y(h_idx,:); Yy = Y(y_idx,:);
th = t(h_idx);   ty = t(y_idx);

ypred = zeros(size(Yy));
tt = [H01_event_PP(32:41).latency]/500;
tt = tt - ty(1); ty = ty - ty(1); th = th - th(1);
figure; hold on;
for tDelta = tt
    tShift = th - tDelta;
    hShift = cell2mat( arrayfun(@(c) ...
            interp1(th, Yh(:,c), tShift, 'nearest', 'extrap'), ...
            1:size(Yh,2), 'UniformOutput',false) );
    plot(th, hShift);
    hDelta = cell2mat( arrayfun(@(c) ...
            interp1(th, hShift(:,c), ty, 'nearest', 'extrap'), ...
            1:size(Yh,2), 'UniformOutput',false) );
    ypred = ypred + hDelta;
end

figure; plot(ypred, Yy, '.');

%% helper functions

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

function T = eventBoundTimes(evs)
    intvl = diff([evs.init_time]); 
    intvlFromPrev = [inf, intvl]; intvlToNext = [intvl, inf];
    %inTrain = (intvlFromPrev < 2) | (intvlToNext < 2);
    %firstOfTrain = (intvl >= 2) & (intvl < 10);
    %prickBefore = intvl >= 10;
    firstOfTrain = (intvlToNext < 2) & (intvlFromPrev >= 2);
    lastOfTrain = (intvlFromPrev < 2) & (intvlToNext >= 2);
    prickBefore = (intvlToNext >= 2) & (intvlToNext < 10) & (intvlFromPrev >= 2);

    evStart = find(firstOfTrain); evEnd = zeros(size(evStart));
    for idx = 1:length(evStart)
        eEnd = find(lastOfTrain);
        eEnd = eEnd(eEnd >= evStart(idx));
        if idx < length(evStart)
            eEnd = eEnd(eEnd <= evStart(idx+1));
        end
        if isempty(eEnd)
            evEnd(idx) = -1;
        else
            evEnd(idx) = max(eEnd);
        end
    end
    evStart = evStart(evEnd >= 0); evEnd = evEnd(evEnd >= 0);

    ev0 = zeros(size(evStart));
    for idx = 1:length(evStart)
        e0 = find(prickBefore);
        e0 = e0(e0 <= evStart(idx));
        if idx > 1
            e0 = e0(e0 >= evStart(idx-1));
        end
        if isempty(e0)
            ev0(idx) = -1;
        else
            ev0(idx) = min(e0);
        end
    end

    T = [ev0; evStart; evEnd]; T = T';
end

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
end