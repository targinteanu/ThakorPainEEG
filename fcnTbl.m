function [tblOut, eegTbl, epochTbl, epochSpecTbl] = ...
    fcnTbl(eegTbl, epochTbl, epochSpecTbl, fcn, vars)
    makeSubtbl = @(tbl, vars) tbl(:, ismember(tbl.Properties.VariableNames, vars));

    if nargin > 4
        eegTbl       = makeSubtbl(eegTbl,       vars);
        epochTbl     = makeSubtbl(epochTbl,     vars);
        epochSpecTbl = makeSubtbl(epochSpecTbl, vars);
    end
    
    tblOut = epochTbl; 
    W = height(epochTbl); H = width(epochTbl); 
    for c = 1:H
        for r = 1:W
            curSpecs = epochSpecTbl{r,c}{:};
            curEpocs = epochTbl{r,c}{:};
            %Y = cell(size(curEpocs)); t = cell(size(Y));
            tY = cell(size(curEpocs));
            for trl = 1:length(curEpocs)
                curSpec = curSpecs{trl};
                curEpoc = curEpocs{trl};
                if ~isempty(curEpoc)
                    [Y_trl,t_trl] = fcn(curSpec, curEpoc);
                    tY{trl} = cat(3,t_trl,Y_trl);
                end
            end
            tblOut{r,c} = tY;
        end
    end
end