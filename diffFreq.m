function [Y, PF] = diffFreq(w,F,bnd, tbl)

    if nargin < 4 
        if ~exist('BandTableHz')
            load('BrainwaveFrequencyTable.mat');
            global BandTableHz
        end
        tbl = BandTableHz;
        if nargin < 3
            bnd = [];
        end
    end

    PF = arrayfun(@(c) peakFreq(w,F(c,:), bnd, tbl), ...
        1:size(F,1));

    Y = abs(PF - PF')./(PF + PF');

end