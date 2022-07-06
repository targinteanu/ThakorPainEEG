function [CoG, maxval, pklocs] = peakFreq(w,P, bnd, tbl)

    if nargin < 4 
        if ~exist('BandTableHz')
            load('BrainwaveFrequencyTable.mat');
            global BandTableHz
        end
        tbl = BandTableHz;
    end

    if nargin < 3
        bnd = [min(w),max(w)]; 
    end
    if isa(bnd,'char') | isa(bnd,'string')
        bnd = band2freqs(bnd, tbl);
    end

    aband = (w >= bnd(1)) & (w <= bnd(2));
    P = P(aband); w = w(aband);
    [~,maxIdx] = max(P); maxval = w(maxIdx);
    [pkP,pklocs] = findpeaks(P, w);
    if length(pklocs) > 1
        % order by magnitude of peak
        [~,ord] = sort(pkP);
        pklocs = pklocs(fliplr(ord));
    end
    CoG = sum(w.*P)/sum(P);
end