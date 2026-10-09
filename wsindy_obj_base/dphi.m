function phifun = dphi(phifun,diff_order)
    syms y; 
    s = diff(phifun(y),diff_order);
%     double(subs(s,y,pi/36))
%     if isAlways(s==1)
%         phifun = @(varargin) varargin{1}*0+1;
%     elseif isAlways(s==0)
%         phifun = @(varargin) varargin{1}*0;
%     else
        phifun = matlabFunction(s);
%     end
end