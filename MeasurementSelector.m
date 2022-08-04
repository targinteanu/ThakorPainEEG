function [fcn, yname, ylims] = MeasurementSelector()

load("BrainwaveFrequencyTable.mat");
global BandTableHz

%% plot type selection 

plotOpts1D =   {'Peak Freq (Hz)', ...
                'Band Mean (\muV^2 s^2)', ...
                'Band Density', ...
                'Node Degree', ...
                'Connectivity Strength', ...
                'Average Frequency Difference', ...
                'Clustering Coefficient', ...
                'Neighbor Average Distance (mm)', ...
                'Neighbors Average ', ...
                'Node Average '};
plotOpts0D =   {'Custom Correlation', ...
                'Frequency Assortativity'};
plotOpts2D =   {'WPLI ', ...
                'PLV ', ...
                'Frequency Difference', ...
                'Channel Distance (mm)', ...
                'Adjacency'};
yname = [];


% basic frequency function definition 
meanAmp = @(w,P,band) mean(P(:, ( (w >= band(1))&(w <= band(2)) )), 2);
ampDensity = @(w,P,band) sum(P(:, ( (w >= band(1))&(w <= band(2)) )), 2) ./ sum(P,2);
basicFreqFcnOpts = {meanAmp, ampDensity};


% initial request for plot type 
plotOpts = [plotOpts1D, plotOpts0D];
[~,PLOTSEL] = listdlg_selectWrapper(plotOpts, 'single', 'Plot What?');


% **** handle correlation ****
isCorrelation = false;
if strcmp(PLOTSEL{1}, 'Custom Correlation')
    yname = 'Correlation Between ';
    isCorrelation = true; 
    PLOTSEL = cell(1,2);
    plotOpts = plotOpts1D;
    for idx = 1:length(PLOTSEL)
        % request each individual variable 
        PS = listdlg_selectWrapper(plotOpts, 'single', ...
            ['Correlation Between: (',num2str(idx),' of ',num2str(length(PLOTSEL)),')']);
        PLOTSEL{idx} = plotOpts{PS};
    end
end


fcn = cell(1,2);
for idx = 1:length(PLOTSEL)
    plotSel = PLOTSEL{idx};

    % request parameters for network variables 
    nwFcn = []; nwSel = '';
    nwFcn_needs_bnd = false;
    Pcut = [];
    if sum(strcmp(plotSel, {...
            'Node Degree', ...
            'Clustering Coefficient', ...
            'Neighbor Average Distance (mm)', ...
            'Frequency Assortativity', ...
            'Adjacency', ...
            'Neighbors Average '}) )
        % percentile cutoff must be selected
        Pcut = inputdlg('Network Cutoff Percentile (%)', ...
            [plotSel,' (',num2str(idx),' of ',num2str(length(PLOTSEL)),') Network Cutoff']);
        Pcut = str2num(Pcut{1});
    end
    if sum(strcmp(plotSel, {...
            'Connectivity Strength', ...
            'Node Degree', ...
            'Clustering Coefficient', ...
            'Neighbor Average Distance (mm)', ...
            'Adjacency', ...
            'Neighbors Average '}) )
        % network method must be specified 
        [~,nwSel] = listdlg_selectWrapper(plotOpts2D, 'single', 'Network Method:');
        nwSel = nwSel{1};
        if strcmp(nwSel, 'WPLI ')
            nwFcn = @(Spec,EEG, bnd) ...
                WPLI(Spec.frequencySpectrum, Pcut, bnd, Spec.frequency2side, BandTableHz);
            nwFcn_needs_bnd = true;
        elseif strcmp(nwSel, 'PLV ')
            nwFcn = @(Spec,EEG, bnd) ...
                PLV();
            nwFcn_needs_bnd = true;
        end
    end

    % *** handle 2D-to-1D functions and neighbor functions ***
    neighborLayers = 0; % number of times to iterate avg_neighbors
    fcn_2D_to_1D = [];  % 2D to 1D conversion function
    fcn_2D_to_1D_needs_nwFcn = false;
    if strcmp(plotSel, 'Neighbors Average ')
        % Neighbors Average Value - of what?
        contChain = true;
        while contChain
            yname = [yname, plotSel];
            plotOpts = [plotOpts1D, plotOpts2D];
            plotSel = listdlg_selectWrapper(plotOpts, 'single', 'Average What?');
            plotSel = plotOpts{plotSel};
            if sum(strcmp(plotSel, plotOpts2D))
                % use avg_neighbor_2D and end the chain
                fcn_2D_to_1D = @(a,b,c,F) avg_neighbor_2D(a,b,c,F);
                fcn_2D_to_1D_needs_nwFcn = true;
                contChain = false;
            else
                % use avg_neighbor and continue the chain only if another average selected
                neighborLayers = neighborLayers + 1;
                contChain = strcmp(plotSel, 'Neighbors Average ');
            end
        end
    end
    if strcmp(plotSel, 'Node Average ')
        % Node Average Value - of what?
        yname = [yname, plotSel];
        fcn_2D_to_1D = @(a,b,c) avg_node_2D(a,b,c);
        plotOpts = plotOpts2D;
        plotSel = listdlg_selectWrapper(plotOpts, 'single', 'Average What?');
        plotSel = plotOpts{plotSel};
    end


    % frequency band must be selected ?
    if (sum(strcmp(plotSel,{...
            'Peak Freq (Hz)', ...
            'Band Mean (\muV^2 s^2)', ...
            'Band Density', ...
            'Frequency Assortativity', ...
            'Average Frequency Difference', ...
            'WPLI ', ...
            'PLV ', ...
            'Frequency Difference'}) ) ) ...
            | nwFcn_needs_bnd
        [bnd,bndname] = pickFrequency([plotSel,' (',num2str(idx),' of ',num2str(length(PLOTSEL)),')']);
        yname = [yname,bndname,' ',plotSel];

        % network 
        if nwFcn_needs_bnd
            nwFcn = @(Spec,EEG) nwFcn(Spec,EEG, bnd);
            nwSel = [bndname,' ',nwSel];
        end
    end
        
    % *** handle 2D-to-1D functions and neighbor functions ***
    if fcn_2D_to_1D_needs_nwFcn
        fcn_2D_to_1D = @(a,b,c) fcn_2D_to_1D(a,b,c, nwFcn);
    end

    isBasicFreqFcn = find(strcmp(plotSel, {...
        'Band Mean (\muV^2 s^2)', ...
        'Band Density'}) );
    if ~isempty(isBasicFreqFcn)
        fcn0 = basicFreqFcnOpts{isBasicFreqFcn};
        if isa(bnd, 'char') | isa(bnd, 'string')
            bnd = band2freqs(bnd, BandTableHz);
        end
        fcn{idx} = @(Spec, EEG) frqFcnEpoch(Spec, EEG, @(w,P) fcn0(w,P,bnd));
        if strcmp(plotSel, 'Band Density')
            ylims = [0 1];
        elseif strcmp(plotSel, 'Band Mean (\muV^2 s^2)')
            ylims = [];
        end

    elseif strcmp(plotSel, 'Peak Freq (Hz)')
        % peak freq
        fcn{idx} = @(Spec, EEG) peakFreqEpoch(Spec, EEG, bnd, BandTableHz);
        if isa(bnd,'char') | isa(bnd,'string')
            ylims = band2freqs(bnd, BandTableHz);
        else
            ylims = bnd;
        end

    elseif strcmp(plotSel, 'Frequency Assortativity')
        % assortativity
        fcn{idx} = @(Spec, EEG) assortativity(Spec, EEG, bnd, Pcut, BandTableHz);
        ylims = [-1, 1];
    elseif strcmp(plotSel, 'Node Degree')
        % node degree
        fcn{idx} = @(Spec, EEG) nodeDegree(Spec, EEG, nwFcn);
        ylims = 'numchan';
    elseif strcmp(plotSel, 'Connectivity Strength')
        % conn strength
        fcn{idx} = @(Spec, EEG) avg_node_2D(Spec, EEG, ...
            nwFcn );
        ylims = [];
    elseif strcmp(plotSel, 'Average Frequency Difference')
        % Freq Diff
        fcn{idx} = @(Spec, EEG) avg_node_2D(Spec, EEG, ...
            @(SO) diffFreq(SO.frequency1side, SO.powerSpectrum, bnd, BandTableHz) );
        ylims = [0, 1];
    elseif strcmp(plotSel, 'Clustering Coefficient')
        % Clustering Coeff
        fcn{idx} = @(Spec, EEG) clusteringCoeff(Spec, EEG, nwFcn);
        ylims = [0, 1];
    elseif strcmp(plotSel, 'Neighbor Average Distance (mm)')
        % Avg Neighbor Dist
        fcn{idx} = @(Spec, EEG) avg_neighbor_2D(Spec, EEG, ...
            @(SO,EO) chanDistance(EO), ...
            nwFcn);
        ylims = [];

    elseif strcmp(plotSel, 'Frequency Difference')
        fcn{idx} = @(Spec, EEG) diffFreq(Spec.frequency1side, Spec.powerSpectrum, bnd, BandTableHz);
        ylims = [0, 1];
    elseif strcmp(plotSel, 'WPLI ')
        fcn{idx} = @(Spec, EEG) WPLI(Spec.frequencySpectrum, Pcut, bnd, Spec.frequency2side, BandTableHz);
        ylims = [0, 1];
    elseif strcmp(plotSel, 'PLV ')
        fcn{idx} = @(Spec, EEG) PLV();
        ylims = [];
    elseif strcmp(plotSel, 'Channel Distance (mm)')
        fcn{idx} = @(Spec, EEG) chanDistance(EEG);
        ylims = [];
    end


    % *** handle 2D-to-1D functions and neighbor functions ***
    if ~isempty(fcn_2D_to_1D)
        fcn{idx} = @(Spec, EEG) fcn_2D_to_1D(Spec, EEG, fcn{idx});
    end
    for l = 1:neighborLayers
        fcn{idx} = @(Spec, EEG) avg_neighbor(Spec, EEG, fcn{idx}, Pcut, bnd, BandTableHz);
    end


    % describe network  
    if ~isempty(Pcut) | ~isempty(nwSel)
        yname = [yname,' in ',nwSel,'Network'];
        if ~isempty(Pcut)
            yname = [yname,' of top ',num2str(100-Pcut),'%'];
        end
    end


    % **** handle correlation ****
    if isCorrelation
        if idx == 1
            yname = [yname,' and '];
        else
            fcn = @(Spec, EEG) fcnCorr(fcn{1}, fcn{2}, Spec, EEG);
            ylims = [-1, 1];
        end
    else
        fcn = fcn{1};
    end
end

%% helper functions 

function [bandrange, bandname] = pickFrequency(ttl)
    freqOpts = [BandTableHz.Properties.RowNames; 'custom'; 'all'];
    bandrange = listdlg_selectWrapper(freqOpts, 'single', [ttl,': Specify Frequency Band']);
    if bandrange == (length(freqOpts) - 1)
        % custom 
        bandrange = inputdlg({'Minimum Frequency (Hz):', 'Maximum Frequency (Hz):'},...
            'Specify Custom Frequency Band:');
        bandname = [bandrange{1},'-',bandrange{2},'Hz'];
        bandrange = arrayfun(@(i) str2double(bandrange{i}), 1:length(bandrange));
    elseif bandrange == length(freqOpts)
        % all frequency 
        bandrange = [];
        bandname = '';
    else
        % named band 
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
            sel = questdlg('Exit?');
            ok = strcmp(sel, 'Yes');
        if ~ok
            [sel, ok] = listdlg('ListString',list, 'SelectionMode',SelectionMode, 'PromptString',PromptString);
        end
    end
    listOut = list(sel);
end

function t = getTimes(epoch_EEG)
    if isempty(epoch_EEG)
        t = [];
    else
        t = arrayfun(@(eeg) mean([eeg.xmin, eeg.xmax]), epoch_EEG);
    end
end

%% time-frequency functions 

function [Y,times] = frqFcnEpoch(epoch_Spec, epoch_EEG, fcn)
    Y = cell2mat( arrayfun(@(eegSpec) ...
        fcn(eegSpec.frequency1side, eegSpec.powerSpectrum), ...
        epoch_Spec, 'UniformOutput', false) );
    times = getTimes(epoch_EEG);
    times = repmat(times,size(Y,1),1);
    times = times'; Y = Y';
end

function [Freq,times] = peakFreqEpoch(epoch_Spec, epoch_EEG, bnd, tbl)
% inspect peakFreq across epochs (time)
    nchan = epoch_Spec(1).nbchan;
    PAFs = zeros(nchan, length(epoch_Spec), 2);
    for chan = 1:nchan
        PAFs(chan,:,1) = arrayfun(@(eegSpec) ...
            peakFreq(eegSpec.frequency1side, eegSpec.powerSpectrum(chan,:), ...
            bnd, tbl), ...
            epoch_Spec);
        PAFs(chan,:,2) = arrayfun(@(eeg) ...
            mean([eeg.xmin, eeg.xmax]), epoch_EEG);
    end
    times = PAFs(:,:,2)'; Freq = PAFs(:,:,1)';
end

function [Y,t] = rawData(EEGObj, chanSelect)
    if nargin < 2
        chanSelect = [];
    end
    if ~isempty(chanSelect)
        EEGObj = pop_select(EEGObj, 'channel', chanSelect);
    end
    Y = EEGObj.data; t = EEGObj.times/1000; 
    t = repmat(t, size(Y,1), 1); 
    t = t'; Y = Y';
end

%% network functions 

function [nd,t] = nodeDegree(SpectObj, EEGObj, networkFcn)
    nd = zeros(SpectObj(1).nbchan,length(SpectObj));
    for s = 1:length(SpectObj)
        [~,A] = networkFcn(SpectObj(s), EEGObj(s));
        nd(:,s) = sum(A);
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(nd,1), 1); 
    nd = nd'; t = t';
end

function [C,t] = clusteringCoeff(SpectObj, EEGObj, networkFcn)
    C = zeros(SpectObj(1).nbchan,length(SpectObj));
    for s = 1:length(SpectObj)
        [~,A] = networkFcn(SpectObj(s), EEGObj(s));
        k = sum(A,2); Ct = zeros(size(k));
        for c = 1:SpectObj(s).nbchan
            Nidx = A(:,c);
            N = A(Nidx, Nidx);
            Ct(c) = sum(N(:));
        end
        C(:,s) = Ct./(k.*(k-1));
        % C(Ct==0,s) = 0; % ??
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(C,1), 1); 
    C = C'; t = t';
end

%% functions for converting dimension or neighbor averaging 

function [A,t] = avg_node_2D(SpectObj, EEGObj, nodeFcn)
    % 2D to 1D
    % each channel's average nodeFcn with all others 
    A = zeros(SpectObj(1).nbchan,length(SpectObj));
    for s = 1:length(SpectObj)
        W = nodeFcn(SpectObj(s), EEGObj(s));
        A(:,s) = mean(W, 'omitnan');
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(A,1), 1); 
    A = A'; t = t';
end

function [A,t] = avg_neighbor_2D(SpectObj, EEGObj, neighborFcn, networkFcn)
    % convert 2D to 1D 
    A = zeros(SpectObj(1).nbchan,length(SpectObj));
    for s = 1:length(SpectObj)
        W = neighborFcn(SpectObj(s), EEGObj(s)); 
        [~,Adj] = networkFcn(SpectObj(s), EEGObj(s));
        for idx = 1:size(Adj,1)
            Adj(idx,idx) = 1;
        end
        A(:,s) = arrayfun(@(c) mean(W(Adj(:,c),c)), 1:size(W,1));
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(A,1), 1); 
    A = A'; t = t';
end

function [A,t] = avg_neighbor(SpectObj, EEGObj, neighborFcn, networkFcn)
    % 1D to 1D
    Y = neighborFcn(SpectObj, EEGObj); Y = Y';
    A = zeros(size(Y));
    for s = 1:length(SpectObj)
        [~,Adj] = networkFcn(SpectObj(s), EEGObj(s));
        A(:,s) = arrayfun(@(c) mean(Y(Adj(:,c),s)), 1:size(Y,1));
    end
    t = getTimes(EEGObj);
    t = repmat(t, size(A,1), 1); 
    A = A'; t = t';
end

%% output a 0D quantity/correlation 

function [Y,t,Y1,Y2] = fcnCorr(fcn1, fcn2, var1, var2)
    [Y1,t1] = fcn1(var1, var2); [Y2,t2] = fcn2(var1, var2);
    Y1 = Y1'; Y2 = Y2'; t1 = t1'; t2 = t2';

    t = sort(unique([t1(:);t2(:)]));
    if length(t) > 1
        Y1 = cell2mat( arrayfun(@(c) ...
            interp1(t1(c,:), Y1(c,:), t, 'linear', 'extrap'), ...
            1:size(Y1,1), 'UniformOutput',false) )';
        Y2 = cell2mat( arrayfun(@(c) ...
            interp1(t2(c,:), Y2(c,:), t, 'linear', 'extrap'), ...
            1:size(Y2,1), 'UniformOutput',false) )';
    end

    Y = zeros(length(t),1);
    for s = 1:length(t)
        Y(s) = corr(Y1(:,s), Y2(:,s), 'Type', 'Spearman', 'Rows', 'complete');
    end
end

function [Af,t] = assortativity(SpectObj, EEGObj, bnd, cutoffPercentile, tbl)
    % 0D
    if nargin < 5
        tbl = BandTableHz;
        if nargin < 4
            cutoffPercentile = [];
            if nargin < 3
                bnd = [];
                if nargin < 2
                    EEGObj = [];
                end
            end
        end
    end
    Af = zeros(size(SpectObj));
    for s = 1:length(SpectObj)
        [~,Adj] = WPLI(SpectObj(s).frequencySpectrum, cutoffPercentile, bnd, ...
            SpectObj(s).frequency2side, tbl);
        [~,PF] = diffFreq(SpectObj(s).frequency1side, SpectObj(s).powerSpectrum, bnd, tbl);
        y = arrayfun(@(c) mean(PF(Adj(:,c))), 1:length(PF));
        Af(s) = corr(PF', y', 'Type','Spearman', 'Rows', 'complete');
    end
    t = getTimes(EEGObj);
    Af = Af'; t = t';
end

%% output a 2D matrix 

function D = chanDistance(EEG)
    chlocs = EEG.chanlocs;
    for c = 1:length(chlocs)
        q = chlocs(c).X + chlocs(c).Y + chlocs(c).Z;
        if isempty(q) | isnan(q)
            chlocs(c).X = nan; 
            chlocs(c).Y = nan; 
            chlocs(c).Z = nan;
        end
    end
    XYZrow = cat(3,[chlocs.X],[chlocs.Y],[chlocs.Z]);
    XYZcol = cat(3,[chlocs.X]',[chlocs.Y]',[chlocs.Z]');
    dXYZ = XYZrow - XYZcol; 
    D = sum(dXYZ.^2, 3).^.5;
end

end