function [phifun,m,p] = get_phi_handle(N,varargin)
% req args: N, opt args: phi_class, tau, tauhat, m, p, k, maxd

    default_phi_class = 1;
    default_tau = 10^-10;
    default_tauhat = 2;
    default_m = floor((N-1)/8);
    default_p = 0;
    default_k = [];
    default_maxd = 4;
    
    inPars = inputParser;
    validScalarPosNum = @(x) isnumeric(x) && all(x >=0);
    addParameter(inPars,'phi_class',default_phi_class);
    addParameter(inPars,'tau',default_tau,validScalarPosNum);
    addParameter(inPars,'tauhat',default_tauhat,validScalarPosNum);
    addParameter(inPars,'m',default_m,validScalarPosNum);
    addParameter(inPars,'p',default_p);
    addParameter(inPars,'k',default_k,validScalarPosNum);
    addParameter(inPars,'maxd',default_maxd,validScalarPosNum);
    parse(inPars,varargin{:});
        
    phi_class = inPars.Results.phi_class;
    tau = inPars.Results.tau;
    tauhat = inPars.Results.tauhat;
    m = inPars.Results.m;
    p = inPars.Results.p;
    k = inPars.Results.k;
    maxd = inPars.Results.maxd;
    
    if isequal(class(phi_class),'function_handle') 
        p = NaN;
        phifun = phi_class;
        if ~isempty(k)
            if isempty(tauhat)
                tauhat = default_tauhat;
            end
            m = get_tf_support(phifun,N,tauhat,k);
        end
    elseif isequal(phi_class,1)
        if ~isempty(k)
            if isempty(tauhat)
                tauhat = default_tauhat;
            end
            if ~(p==0)
%                 l = @(m,k,N) log((2*m-1)./m.^2).*((4*pi^2*k^2*m.^2-3*N^2*tauhat^2)-2*N^2*tauhat^2*p);
%                 mstar1 = sqrt(3)/pi*N/2/k*tauhat;
%                 mstar2 = 1/pi*tauhat*(N/2)/k*sqrt(3+8*p);
%                 m = floor(min(fzero(@(m)l(m,k,N), [mstar1 mstar2]),(N-1)/2));
%               
                phifun = @(x) (1-x.^2).^p;
                m = get_tf_support(phifun,N,tauhat,k);

                tau = ((2*m-1)./m.^2)^p;
            else
                if isempty(tau)
                    tau = default_tau;
                end
                l = @(m,k,N) log((2*m-1)./m.^2).*(4*pi^2*k^2*m.^2-3*N^2*tauhat^2)-2*N^2*tauhat^2*log(tau);
                mstar1 = sqrt(3)/pi*N/2/k*tauhat;
                mstar2 = 1/pi*tauhat*(N/2)/k*sqrt(log(exp(1)^3/tau^8));
                m = floor(min(fzero(@(m)l(m,k,N), [mstar1 mstar2]),(N-1)/2));
                p = max(maxd+1,ceil(log(tau)/log(1-(1-1/m)^2)));
            end
        elseif isempty(p)
            p = ceil(max(log(tau)/log((2*m-1)/m^2),maxd+1));
        else 
        end
        phifun = @(x) (1-x.^2).^p;
    elseif isequal(phi_class,2)
        if ~isempty(k)
            m = floor(min(1+N*tauhat/2/pi/k*sqrt(-2*log(tau)),(N-1)/2));
            p = 2*pi*k/tauhat/N;
        elseif ~isempty(tau)
            a = sqrt(-2*log(tau));
            p = (1-1/m)/a;
        end
        phifun = @(x) exp(-(m*p*x).^2/2); % p = dx/sig
    end 
end