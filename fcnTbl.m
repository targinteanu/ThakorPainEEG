function [tblOut, tblIn] = ...
    fcnTbl(tblIn, fcn, vars)
    makeSubtbl = @(tbl, vars) tbl(:, ismember(tbl.Properties.VariableNames, vars));

    if nargin > 4
        tblIn = makeSubtbl(tblIn, vars);
    end
    
    tblSpec = tblIn(4,:);
    tblEpoc = tblIn(3,:);
    tblOut = tblEpoc; tblOut.Properties.RowNames = {'tY'};
    W = width(tblIn); 
    for c = 1:W
        curSpecs = tblSpec{1,c}{:};
        curEpocs = tblEpoc{1,c}{:};
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
        tblOut{1,c} = {tY};
    end

    tblOut = [tblIn; tblOut];
end