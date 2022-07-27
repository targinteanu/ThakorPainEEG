function [W, A] = WPLI(Z, cutoffPercentile, bnd, w, tbl)

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
        C(r1,:,r2) = Z(r1,:).*conj(Z(r2,:));
    end
end

C = imag(C);
W = abs(mean(C,2))./mean(abs(C),2);
W = squeeze(W);

Wval = sort(W(:)); cutoff = Wval(round((cutoffPercentile/100)*length(Wval)));
if ~sum(Wval < cutoff)
    % exclude ties, or else every node will be paired 
    A = W > cutoff;
else
    % allow ties to avoid situation in which no nodes are paired 
    A = W >= cutoff; 
end

end