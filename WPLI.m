function [W, A] = WPLI(Z, cutoffPercentile)

if nargin < 2
    cutoffPercentile = [];
end
if isempty(cutoffPercentile)
    cutoffPercentile = 60;
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
A = W >= cutoff; 

end