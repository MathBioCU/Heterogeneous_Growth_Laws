%%% even case - cosine basis
function [f,fp] = optTFcos(p,d)
    k=p+1+d;
    r_center = ones(1,k+1);
    r_derivs = [];
    for r=0:2:2*p
        r_derivs =[r_derivs; arrayfun(@(i) (-1)^i*(pi*i)^r*(-mod(r+1,2))^(ceil(r/2)),0:k)];
    end
    V = [r_center;r_derivs];
    c = lsqminnorm(V,[1;zeros(size(V,1)-1,1)]);
    f = @(t) cospolyval(c,pi*(0:k),t);
    fp = @(t) -sinpolyval(pi*(0:k)'.*c,pi*(0:k),t);
end

function y = cospolyval(c,fs,x)
    y = x*0;
    for i=1:length(fs)
        y = y + c(i)*cos(fs(i)*x);
    end
end

function y = sinpolyval(c,fs,x)
    y = x*0;
    for i=1:length(fs)
        y = y + c(i)*sin(fs(i)*x);
    end
end
