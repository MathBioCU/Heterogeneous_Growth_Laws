function x = ismember_rows(a,b)
    if isempty(a)
        x = [];
    elseif isempty(b)
        x = false(size(a));
    else
        x = ismember(a,b,'rows');
    end
end