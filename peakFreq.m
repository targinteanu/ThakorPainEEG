function [CoG, maxval, pklocs] = peakFreq(w,P, bnd)
    if nargin < 3
        bnd = [min(w),max(w)]; 
    end
    if isa(bnd,'char') | isa(bnd,'string')
        if strcmpi(bnd,'alpha') | strcmpi(bnd,'a')
            bnd = [9,11];
        elseif strcmpi(bnd,'beta') | strcmpi(bnd,'b')
            bnd = [13,30];
        elseif strcmpi(bnd,'theta') | strcmpi(bnd,'t')
            bnd = [4,8];
        elseif strcmpi(bnd,'gamma') | strcmpi(bnd,'g')
            bnd = [30,80];
        elseif strcmpi(bnd,'delta') | strcmpi(bnd,'d')
            bnd = [.5,4];
        else
            bnd = [];
        end
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