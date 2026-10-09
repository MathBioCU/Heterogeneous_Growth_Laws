

function [mts,pts,sig_ests,corner] = findcorners(xobs,t,tau,tauhat,phi_class)
    T = length(t);
    
    [corner,sig_ests] = findcornerpts(xobs,t);
    k = corner(2);
    if isequal(phi_class,1)
        l = @(m,k,N) log((2*m-1)./m.^2).*(4*pi^2*k^2*m.^2-3*N^2*tauhat^2)-2*N^2*tauhat^2*log(tau);
        mnew = fzero(@(m)l(m,k,T), [1 2/sqrt(tau)]);
        if mnew>T/2-1
            mnew = T/2/k;
        end
        mts = min(floor((T-1)/2),ceil(mnew)); 
        pts= max(2,floor(log(tau)/log(1-(1-1/mts)^2)));
    elseif isequal(phi_class,2)
        mnew = 1+T*tauhat/2/pi/k*sqrt(-2*log(tau));
        mts = min(floor((T-1)/2),ceil(mnew)); 
        pts= 2*pi*k/tauhat/T;
    elseif isequal(class(phi_class), 'function_handle')
        mts = get_tf_support(phi_class,T,tauhat,k);
        pts = NaN;
    end
    
end