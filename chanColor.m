function rgb = chanColor(chloc, chlocs)
% Assign a unique color to channel based on location. Channels will form a
% color gradient depending on their X, Y, Z coordinates. 
% 
% Inputs: 
%   chloc: object of current channel to assign color 
%   chlocs: array of other channel objects being studied 
% Outputs: 
%   rgb: [r, g, b] color triplet 

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
    rgb = .7 * rgb;
end