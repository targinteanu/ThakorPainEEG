function [freqs,bandname] = band2freqs(bandname, tbl)

    if nargin < 2 
        if ~exist('BandTableHz')
            load('BrainwaveFrequencyTable.mat');
            global BandTableHz
        end
        tbl = BandTableHz;
    end
    
    try 
        freqs = [tbl.minimumFrequency(bandname), tbl.maximumFrequency(bandname)];
    catch
        RowNames = tbl.Properties.RowNames;
        tf = strcmpi(bandname,RowNames);
        if ~(sum(tf) == 1)
            firstLetter = arrayfun(@(i) RowNames{i}(1), 1:length(RowNames));
            tf = bandname(1) == firstLetter;
            if ~(sum(tf) == 1)
                tf = lower(bandname(1)) == firstLetter;
                if ~(sum(tf) == 1)
                    tf = upper(bandname(1)) == firstLetter;
                end
            end
        end
        if (sum(tf) == 1)
            bandname = find(tf);
            freqs = [tbl.minimumFrequency(bandname), tbl.maximumFrequency(bandname)];
        else
            % error(['Unrecognized frequency band ',bandname]);
            freqs = [];
        end
    end

end

