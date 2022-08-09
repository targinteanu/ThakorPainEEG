%% Start eeglab
clear
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Set Filepaths
home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)

%%
fp = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG/Data_Chronic Pain/Preprocessed 2022-08-06 19.51.25/Postprocessed 2022-08-07 01.56.52 -- 4s epoch, 3.8s overlap/';
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

%%
clear H01 P04

%% select what to plot
[fcn, yname, ylims] = MeasurementSelector();

%%
[Y,t] = fcn(H01_Epoch_w_PP, H01_Epoch_t_PP); t = t(:,1);
BL = fcn(H01_w_BL, H01_t_BL); 
BLe = fcn(H01_Epoch_w_BL, H01_Epoch_t_BL); 
% Y = Y - BL; % remove baseline 
Y = (Y - BL)./(std(BLe)/sqrt(size(BLe,1))); % t statistic 
chloc = H01_t_BL.chanlocs;
%Y = Y(:,14); chloc = chloc(14); % cz

%{
flt = designfilt('highpassiir', 'SampleRate', 1/mean(diff(t)), ...
    'PassbandFrequency', .25, 'StopbandFrequency', .15);
Y = filtfilt(flt, Y);
%}
%%
trl = 2;
tsplit = H01_boundTimes_PP(trl,:);
tt = [H01_event_PP(tsplit).latency]/500;
hbnd = tt(1:2) + [0,-.5]; ybnd = tt(2:3);
h_idx = (t <= hbnd(2)) & (t >= hbnd(1));
y_idx = (t <= ybnd(2)) & (t >= ybnd(1));
Yh = Y(h_idx,:); Yy = Y(y_idx,:);
th = t(h_idx);   ty = t(y_idx);

figure; plot(th, Yh); xlabel('time (s)'); ylabel(yname);
title(['Subject H01 PinPrick Trial ',num2str(trl),' - Impulse Response']);

ypred = zeros(size(Yy));
tt = [H01_event_PP(tsplit(2):tsplit(3)).latency]/500;
tt = tt - ty(1); ty = ty - ty(1); th = th - th(1);
%figure; hold on;
for tDelta = tt
    tShift = ty - tDelta;
    hShift = cell2mat( arrayfun(@(c) ...
            interp1(th, Yh(:,c), tShift, 'nearest', 'extrap'), ...
            1:size(Yh,2), 'UniformOutput',false) );
    %plot(th, hShift);
    %{
    hDelta = cell2mat( arrayfun(@(c) ...
            interp1(th, hShift(:,c), ty, 'nearest', 'extrap'), ...
            1:size(Yh,2), 'UniformOutput',false) );
    %}
    %plot(ty, hShift);
    %plot(tDelta, 0, 'vr', 'LineWidth', 2);
    ypred = ypred + hShift;
end
%plot(ty, Yy);
%{
figure; subplot(2,1,1); plot(ypred, Yy, '.');
subplot(2,1,2); plot(ty, Yy); hold on; plot(ty, ypred);
%}

%%
sttl = ['Subject H01 PinPrick Trial ',num2str(trl),' - ',yname];
fig(2) = figure; sgtitle(sttl);
fig(1) = figure; sgtitle(sttl);
H = floor(sqrt(size(Yy,2)));
W = ceil(size(Yy,2)/H);
RP = zeros(2, size(Yy,2));
figure(fig(1));
for idx = 1:size(Yy,2)
    [RP(1,idx), RP(2,idx)] = corr(ypred(:,idx), Yy(:,idx));
    subplot(H,W,idx); 
    plot(ypred(:,idx), Yy(:,idx), '.');
    xlabel('LTI prediction'); ylabel('actual'); 
    title([chloc(idx).labels,...
        ': \rho = ',num2str(RP(1,idx),2),...
        '; p = ',num2str(RP(2,idx),1)]);
end
figure(fig(2));
for idx = 1:size(Yy,2)
    subplot(H,W,idx); 
    plot(ty, ypred(:,idx)); hold on; plot(ty, Yy(:,idx));
    xlabel('time (s)') 
    title([chloc(idx).labels,...
        ': \rho = ',num2str(RP(1,idx),2),...
        '; p = ',num2str(RP(2,idx),1)]);
end
legend({'LTI prediction', 'actual'});

%%
figure; sgtitle(sttl); 
subplot(1,2,1); 
topoplot(RP(1,:), chloc, 'maplimits',[-1 1], 'electrodes','labels');
colorbar;
title('\rho');
subplot(1,2,2); 
topoplot(RP(2,:), chloc, 'maplimits',[0 1],  'electrodes','labels');
colorbar;
title('p');

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