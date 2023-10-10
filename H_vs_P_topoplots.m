% rename to H_vs_P_vs_Baseline_topoplots ?
% create new function: H_vs_P_2x2_topoplots ?

%% Start eeglab
clear
%eeglabpath = 'C:\Program Files\MATLAB\R2022a\eeglab2023.0';
eeglabpath = '/Applications/MATLAB_R2021b.app/toolbox/eeglab2022.0';
addpath(eeglabpath)
eeglab

%% Set Filepaths
home = '/Users/torenarginteanu/Documents/MATLAB/ThakorPainEEG';
cd(home)
load("BrainwaveFrequencyTable.mat");
global BandTableHz
postproDir = uigetdir;

cd(postproDir);
scanfiles = dir; 
scanfiles = scanfiles((~strcmp({scanfiles.name},'.')) & ...
                      (~strcmp({scanfiles.name},'..')) & ...
                      (~strcmp({scanfiles.name},'.DS_Store')) );
scanfiles = scanfiles(~[scanfiles.isdir]);
cd(home); addpath(postproDir);
[~,scanfilesH] = listdlg_selectWrapper({scanfiles.name},'multiple','Select Controls');
[~,scanfilesP] = listdlg_selectWrapper({scanfiles.name},'multiple','Select Patients');
scanfiles = {scanfilesH, scanfilesP};
scanfileNames = {'Control', 'Patient'};
maxNgrp = max( arrayfun(@(s) length(scanfiles{s}), 1:length(scanfiles)) );
maxNgrp = maxNgrp + 1; % make room for combo subj

svloc = [postproDir,'/Comparison Topoplots ',...
    datestr(datetime, 'yyyy-mm-dd HH.MM.SS')];
mkdir(svloc); svloc = [svloc,'/'];

%% select what to plot

nMeas = inputdlg('How many measurements?');
nMeas = str2num(nMeas{1});

fcns = cell(1, nMeas); 
ynames = fcns; ylims = fcns;
for n = 1:nMeas
    [fcns{n}, ynames{n}, ylims{n}] = MeasurementSelector();
end
clear n

%% calculations 
testVars = {'BaselieOpenBefore', 'TempStim', 'PinPrick'};
DATATABLES = cell(size(scanfiles));

for s = 1:length(scanfiles)
    sf = scanfiles{s};

    dataTables = cell(length(sf),nMeas);  
    for subj = 1:length(sf)
        fn = sf{subj}
        load(fn);

        for n = 1:nMeas
        fcn = fcns{n};

        % run calculations on desired variables
        disp('calculating desired variables')
        % construct data tables *************
        
        
        tY_table = fcnTbl(Combined_table, fcn, testVars);

        dataTables{subj, n} = tY_table;
        clear tY_table

        clear fcn
        end
        clear EEG_table Epoch_table EpochSpec_table Spec_table

    end

    DATATABLES{s} = dataTables;
    clear sf dataTables
end
clear testVars

%% channel selection 
%%{
[chansel, chanselName, allchan0, allchan] = ChannelSelector(scanfiles, DATATABLES);
%}
%% combination "subjects" 
varnames = DATATABLES{1}{1,1}.Properties.VariableNames;
comboSubjTbls = cell(1, nMeas);

somechan = true(size(allchan));
for s = 1:length(scanfiles)
    dataTables = DATATABLES{s};
    for subj = 1:size(dataTables,1)
        dataTable = dataTables{subj,1}; % Combined_table with tY
        for c = 1:width(dataTable)
            EEG_all = dataTable{1,c}{1}; % EEG 
            if ~isempty(EEG_all)
                EEG_all = EEG_all(1);
                for ch = 1:length(somechan)
                    somechan(ch) = somechan(ch) & sum( ...
                        strcmpi(allchan(ch).labels, {EEG_all.chanlocs.labels}) );
                end
            end
        end
    end
end
somechan = allchan(somechan);
% somechan = only those selected AND in all EEGs of all subjects of all groups

for n = 1:nMeas
comboSubjTbl = table('size',[length(scanfileNames),length(varnames)], ...
                     'VariableTypes',repmat("cell",size(varnames)), ...
                     'VariableNames',varnames,'RowNames',scanfileNames);

for s = 1:length(scanfiles)

    dataTables = DATATABLES{s};
    cumuTYs = {};

    for subj = 1:size(dataTables,1)

        dataTable = dataTables{subj,n};
        cumuTYsubj = cell(1, width(dataTable));
        for c = 1:width(dataTable)
            EEG_all  = dataTable{1,c}{1};
            if ~isempty(EEG_all)
                EEG_all = EEG_all(1);
                tY_by_trial = dataTable{5,c}{1};
                cumuTY = [];
                for trl = 1:length(tY_by_trial)
                    tY = tY_by_trial{trl};
                    ord = zeros(1,length(somechan));
                    for ch = 1:length(somechan)
                        ord(ch) = find( ...
                            strcmpi(somechan(ch).labels, {EEG_all.chanlocs.labels}) );
                    end
                    cumuTY = cat(1, cumuTY, tY(:,ord,:)); % all trials concatenated
                end
                cumuTYsubj{c} = cumuTY; % one subj; c = stim type
            end
        end
        cumuTYs = [cumuTYs; cumuTYsubj]; % #subjs x #stimtypes

    end

    % concat all subjs (collapse cumuTYs vertically)
    for c = 1:size(cumuTYs,2)
        cumuTY = [];
        for subj = 1:size(cumuTYs,1)
            cumuTY = cat(1, cumuTY, cumuTYs{subj, c});
        end
        comboSubjTbl{s,c} = {cumuTY};
    end

    clear cumuTY cumuTYsubj dataTable ...
          EEG_all tY_by_trial tY ord

    clear dataTables cumuTYs
end
comboSubjTbls{n} = comboSubjTbl;
% everything in comboSubj should be ordered according to somechan
clear comboSubjTbl
end

%% comparison 
% make more robust with more than 2 experimental groups / more baselines? 
varnames = DATATABLES{1}{1}.Properties.VariableNames;
p_alpha = 0.01; % uncertainty for stat significance 
maxstatval = -Inf; minstatval = Inf;
statsTables = cell(size(comboSubjTbls));

for n = 1:nMeas
comboSubjTbl = comboSubjTbls{n};

statsTable = table('size',[3,length(varnames)], ...
                   'VariableTypes',repmat("cell",size(varnames)), ...
                   'RowNames',{'P vs H','H vs baseline','P vs baseline'}, ...
                   'VariableNames',varnames);

% statistical testing 
for c = 1:width(comboSubjTbl)

    s = repmat(struct('chan',somechan(1), 'tstat',0, 'pval',1), ...
               [length(somechan),height(statsTable)]);
    for ch = 1:length(somechan)
        % P vs H
        if (~isempty(comboSubjTbl{1,c}{1}))&(~isempty(comboSubjTbl{2,c}{1}))
            [~,p,~,S] = ...
                ttest2( comboSubjTbl{2,c}{1}(:,ch,2), ...
                comboSubjTbl{1,c}{1}(:,ch,2), ...
                'Vartype', 'unequal' );
            s(ch,1).tstat = S.tstat;
            s(ch,1).pval  = p;
        end
        s(ch,1).chan  = somechan(ch);

        % H vs baseline 
        if ~isempty(comboSubjTbl{1,c}{1})
            [~,p,~,S] = ...
                ttest2( comboSubjTbl{1,c}{1}(:,ch,2), ...
                comboSubjTbl{1,1}{1}(:,ch,2), ...
                'Vartype', 'unequal' );
            s(ch,2).tstat = S.tstat;
            s(ch,2).pval  = p;
        end
        s(ch,2).chan  = somechan(ch);

        % P vs baseline 
        if ~isempty(comboSubjTbl{2,c}{1})
            [~,p,~,S] = ...
                ttest2( comboSubjTbl{2,c}{1}(:,ch,2), ...
                comboSubjTbl{2,1}{1}(:,ch,2), ...
                'Vartype', 'unequal' );
            s(ch,3).tstat = S.tstat;
            s(ch,3).pval  = p;
        end
        s(ch,3).chan  = somechan(ch);

        clear p S
    end

    statsTable{1,c} = {s(:,1)};
    statsTable{2,c} = {s(:,2)};
    statsTable{3,c} = {s(:,3)};
    clear s
end

% for plotting: find bounds and significance 
for r = 1:height(statsTable)
    for c = 1:width(statsTable)
        S = statsTable{r,c}{1};
        %{
        for ch = 1:length(S)
            if S(ch).pval < p_alpha
                S(ch).chan.labels = '*';
            else
                S(ch).chan.labels = '.';
            end
        end
        %}
        statsTable{r,c} = {S};
        S = [S.tstat];
        maxstatval = max(maxstatval, max(S)); minstatval = min(minstatval, min(S));
    end
end

% rename table for plotting 
statsTable.Properties.VariableNames = {...
    'Baseline',          ... BaselineOpenBefore
    'After',             ... BaselineOpenAfter
    'Eyes Closed',       ... BaselineClosedBefore
    'Eyes CLosed After', ... BaselineClosedAfter
    'Heat Stimulus',     ... TempStim
    'Sharp Stimulus',    ... PinPrick 
    'Pressure Stimulus', ... Pressure
    'Sharp + CPM',       ... PinPrickCPM
    'Pressure + CPM'     ... PressureCPM
    };
statsTables{n} = statsTable;

clear S statsTable
end

maxstatval = max(abs(maxstatval), abs(minstatval));
minstatval = -maxstatval;

%% plot comparison 
plot_vs_baseline = {'Heat Stimulus'};
plot_P_vs_H      = {'Baseline'};
flipRC = true;

W = 2*length(plot_vs_baseline) + length(plot_P_vs_H);
fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]); 
%sgtitle([' t statistic; * p < ',num2str(p_alpha)]);

plot_P_vs_H_2 = plot_P_vs_H;
for idx = 1:length(plot_P_vs_H_2)
    plot_P_vs_H_2{idx} = [plot_P_vs_H_2{idx},' Patient vs Control'];
end
plot_P_vs_H_2 = plot_P_vs_H_2';
plot_vs_baseline_2 = cell(1, length(plot_vs_baseline), length(scanfileNames));
for c = 1:length(plot_vs_baseline)
    for r = 1:length(scanfileNames)
        plot_vs_baseline_2{r,c} = [scanfileNames{r},' ',plot_vs_baseline{c},' vs Baseline'];
    end
end
plot_vs_baseline_2 = plot_vs_baseline_2(:);

lblSubplots(nMeas,W, ynames,[plot_P_vs_H_2; plot_vs_baseline_2], ...
            flipRC, 1, 'westoutside');

for n = 1:nMeas
statsTable = statsTables{n};
w = 1;

for idx = 1:length(plot_P_vs_H)
    subplot_Wrapper(nMeas,W, n,w, flipRC);
    pltname = plot_P_vs_H{idx};
    S = statsTable(1, strcmp(pltname, statsTable.Properties.VariableNames));
%    title([pltname,' ',S.Properties.RowNames{1}]);
    S = S{1,1}{1};
    topoplot([S.tstat], [S.chan], ...
        'emarker2', {find([S.pval]<p_alpha), '*', 'k'}, ...
        ...'pmask', [S.pval]<p_alpha, ...
        ...'numcontour',[p_alpha p_alpha], 'contourvals',[S.pval], ...
        'maplimits', [minstatval, maxstatval]); 
    w = w + 1;
end
for idx = 1:length(plot_vs_baseline)
    pltname = plot_vs_baseline{idx};
    S = statsTable([2,3], strcmp(pltname, statsTable.Properties.VariableNames));
    for idx2 = 1:height(S)
        subplot_Wrapper(nMeas,W, n,w, flipRC);
%        title([pltname,' ',S.Properties.RowNames{idx2}]);
        SS = S{idx2,1}{1};
        topoplot([SS.tstat], [SS.chan], ...
            'emarker2', {find([SS.pval]<p_alpha), '*', 'k'}, ...
            'maplimits', [minstatval, maxstatval]); 
        w = w + 1;
    end
end
%colorbar('Location', 'eastoutside');

clear S SS statsTable
end

%saveas(fig, [svloc,yname,' StatsTopoPlot'], 'fig');

%% helper functions 

function subtbl = makeSubtbl(tbl, vars)
    subtbl = tbl(:, ismember(tbl.Properties.VariableNames, vars));
end

function subplot_Wrapper(H,W,r,c,flipRC)
    if flipRC
        HH = H; H = W; W = HH;
        rr = r; r = c; c = rr;
    end
    W = W+1; c = c+1; % pad an extra column for vertical labels 
    p = (r-1)*W + c;
    subplot(H,W,p);
end
function lblSubplots(H,W,rttl,cttl,flipRC,cbCol,cbLcn)
    if nargin < 7
        cbLcn = 'westoutside';
        if nargin < 6
            cbCol = 1;
        end
    end
    %clbl = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

    % label all plots 
    for r = 1:H
        for c = 1:W
            subplot_Wrapper(H,W,r,c,flipRC)
            title([clbl(c),'.',rlbl(r)]);
        end
    end

    if flipRC 
        for r = 1:W
            % label rows
            subplot_Wrapper(W,H,r,0,false)
            ylabel([clbl(r),') ',cttl{r}]);
            % colorbars 
            subplot_Wrapper(W,H,r,cbCol,false)
            colorbar('Location',cbLcn);
        end
        for c = 1:H
            % label cols 
            subplot_Wrapper(W,H,1,c,false)
            title({[rlbl(c),') ',rttl{c}], ...
                [clbl(1),'.',rlbl(c)]})
        end
    else 
        for r = 1:H
            % label rows
            subplot_Wrapper(H,W,r,0,false)
            ylabel([rlbl(r),') ',rttl{r}]);
            % colorbars 
            subplot_Wrapper(H,W,r,cbCol,false)
            colorbar('Location',cbLcn);
        end
        for c = 1:W
            % label cols 
            subplot_Wrapper(H,W,1,c,false)
            title({[clbl(c),') ',cttl{c}], ...
                [clbl(c),'.',rlbl(1)]})
        end
    end

    function let = clbl(numval)
        numval = numval - 1;
        base = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
        protoLet = dec2base(numval, length(base)); let = protoLet;
        shiftbase = ['0123456789',base];
        for idx = 1:length(protoLet)
            loc = find(shiftbase == protoLet(idx));
            let(idx) = base(loc);
        end
    end
    function num = rlbl(numval)
        base = 'ivxlcdm';
        baseval = [1, 5, 10, 50, 100, 500, 1000];
        almostval = baseval - baseval';

        num = '';
        while numval > 0
            % "irregular" cases (iv, ix, etc)
            for idxNext = fliplr(1:length(baseval))
                for idxDiff = 1:length(baseval)
                    if numval == almostval(idxDiff, idxNext)
                        num = [num,base(idxDiff),base(idxNext)];
                        numval = 0;
                        break; 
                    end
                end
                if numval == 0
                    break;
                end
            end
            if numval == 0
                break;
            end

            % "regular" cases (no iv, ix, etc)
            idxFirst = find(numval >= baseval); idxFirst = idxFirst(end);
            num = [num,base(idxFirst)]; numval = numval - baseval(idxFirst);
        end
    end
end