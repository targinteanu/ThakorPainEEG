function outArr = UnwrapCat(Arr1, Arr2)
% unwrap two cell arrays and concatenate them together 
% TO DO: can this be included as a helper of tblReorg? Used elsewhere?
    Arr1 = Arr1{:}; Arr2 = Arr2{:};
    if isempty(Arr1)
        outArr = Arr2;
    elseif isempty(Arr2)
        outArr = Arr1;
    else
        outArr = [Arr1, Arr2];
    end
end