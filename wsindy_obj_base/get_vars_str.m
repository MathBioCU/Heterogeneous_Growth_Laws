function [vars,expr] = get_vars_str(s)
    a = strfind(s,'(');
    b = strfind(s,')');
    expr = s(1:b(1));
    vars = strsplit(s(a(1)+1:b(1)-1),',');
end
