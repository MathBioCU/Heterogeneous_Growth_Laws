
function s=fcn_rep_vars(s,vn)
    [vars,expr] = get_vars_str(s);
    s = strrep(s,expr,'');
    for j=1:length(vn)
        s = strrep(s,vars{j},vn{j});
    end
end