% p = 3; d =0;
% f = optTFpoly(p,d);
% f(linspace(-1,1,100)')

function f = optTFpoly(p,d)
    k=p+1+d;
    r_center = [1 zeros(1,k)];
    r_derivs = zeros(p+1,k+1);
    for r=0:p
        r_derivs(r+1,:) =arrayfun(@(i) prod(i:-1:(i-r+1)),0:2:2*k);
    end
    V = [r_center;r_derivs];
    c = lsqminnorm(V,[1;zeros(p+1,1)]);
    cp = [1;0]*c'; cp = cp(:);
    f = @(t) horn(flipud(cp),t);
end

function y=horn(c,x)
    l = length(c);
    y = c(1);
    for j=1:l-1
        y = y.*x+c(j+1);
    end
end