function [P, A] = PLV(Z, cutoffPercentile, bnd, w, tbl)

if nargin < 5
    tbl = [];
    if nargin < 3
        bnd = [];
        if nargin < 2
            cutoffPercentile = [];
        end
    end
end
if isempty(cutoffPercentile)
    cutoffPercentile = 60;
end

if ~isempty(bnd)
    if isempty(tbl)
        if ~exist('BandTableHz')
            load('BrainwaveFrequencyTable.mat');
            global BandTableHz
        end
        tbl = BandTableHz;
    end
    if isa(bnd,'char') | isa(bnd,'string')
        bnd = band2freqs(bnd, tbl);
    end

    inband = (w >= bnd(1)) & (w <= bnd(2));
    Z = Z(:, inband); 
end

C = zeros([size(Z), size(Z,1)]); 
for r1 = 1:size(Z,1)
    for r2 = 1:size(Z,1)
        C(r1,:,r2) = (Z(r1,:)./abs(Z(r1,:))) ./ (Z(r2,:)./abs(Z(r2,:)));
    end
end

P = abs(mean(C,2)); P = squeeze(P);

Pval = sort(P(:)); cutoff = Pval(round((cutoffPercentile/100)*length(Pval)));
if ~sum(Pval < cutoff)
    % exclude ties, or else every node will be paired 
    A = P > cutoff;
else
    % allow ties to avoid situation in which no nodes are paired 
    A = P >= cutoff; 
end

end