%% choose spatial grid
M = 1000; % number of observed data points
x = linspace(0,2.5*pi,M);

%% choose function to differentiate
k = 5; % number of sine features
max_freq = max(k,ceil(1/32*M/2)); % maximum mode as fraction of nyquist 
c = randn(1,k); % random sine coefficients
fq = randperm(max_freq,k); 
f = @(t) sinpolyval(c,fq,t);
fp = @(t) cospolyval(fq.*c,fq,t);

%% add noise
sig = 0.05; % noise level
Y_cl = f(x');
Y = Y_cl+std(Y_cl)*sig*randn(size(Y_cl));
Yp = fp(x');
U = wsindy_data(Y,x);

figure(10)
plot(x,f(x),'-',x,Y,'ro')
legend({'true fcn','data'})

%% toggle methods to display

toggle_gb = 1; % run global basis differentiation
toggle_gp = 1; % run Gaussian process differentiation
toggle_fd = 1; % run finite difference differentiation
toggle_lp = 1; % run local basis differentiation
toggle_wf = 1; % run weak-form differentiation

%% FD
if toggle_fd==1
    tic;
    fd_ord = 4; % must be even
    fd_Op = get_loc(fd_ord,x,'basis','poly');
    Yp_fd = fd_Op*Y;
    fd_time = toc;
    figure(1); clf
    semilogy(abs(Yp-Yp_fd)/mean(abs(Yp)))
    title('Finite difference')
else
    Yp_fd = Yp;
end

%% glob proj
if toggle_gb==1
    tic;
    gb_ord = ceil(1/4*M/2); %%% set to M/2 + trig basis to Fourier differentiation
    [Globalderiv,Proj_mat,Phi_mat] = get_glob(gb_ord,x,'basis','trig');
    Y_glob = Proj_mat*Y; %%% can also just use as a smoother
    Yp_glob = Globalderiv*Y;
    glob_time = toc;
    figure(2); clf
    semilogy(abs(Yp-Yp_glob)/mean(abs(Yp)))
    title('Global projection')
else
    Yp_glob = Yp;
    Y_glob = Y;
end

%% loc proj
if toggle_lp==1
    tic;
    lp_deg = 3; % local poly basis degree
    U.getcorners; % get FFT corner estimate
    nloc = 4*ceil(lp_deg/2); % half-width of lp stencil
    xloc = x(1:2*nloc+1); % local stencil
    [Mat,PMat] = get_loc_lsq(lp_deg,xloc,M,'basis','poly');
    Ylp_loc = Mat*Y; % deriv estimation
    Yl_loc = PMat*Y; % state estimation
    ploc_time = toc;
    figure(3); clf
    semilogy(abs(Yp-Ylp_loc)/mean(abs(Yp)))
    title('local projection')
else
    Ylp_loc = Yp;
    Yl_loc = Y;
end

%% GP 
if all([toggle_gp==1 toggle_gb==1 toggle_fd==1])
    tic;
    K = (Phi_mat*Phi_mat'); %%% requires glob_proj
    sigest = U.estimate_sigma{1}; %%% estimate noise in data
    Yp_gp = K*((K+sigest^2*(fd_Op'*fd_Op)) \ Yp_fd); %%% requires fd
    gp_time = fd_time+toc;
    figure(4); clf
    semilogy(abs(Yp-Yp_gp)/mean(abs(Yp)))
    title('Gaussian process projection')
else
    Yp_gp = Yp;
end

%% weak deriv
if toggle_wf==1
    tic,
    Dwf = diffOp(1,'meth','wffd');
    U = wsindy_data(Y,x);
    Yp_wf = Dwf.evalterm(U);
    wf_time = toc+fd_time;
    figure(5); clf
    semilogy(abs(Yp-Yp_wf)/mean(abs(Yp)))
    title('Weak derivative')
else
    Yp_wf = Yp;
end

%% diplay results

timez = {fd_time,glob_time,ploc_time,gp_time,wf_time};
errz = num2cell([vecnorm(Yp - [Yp_fd Yp_glob Ylp_loc Yp_gp Yp_wf])]/norm(Yp));
fprintf('runtimes:\n (FD) %.2e;\n (GB) %.2e;\n (LB) %.2e;\n (GP) %.2e;\n (WF) %.2e;\n',timez{:})
fprintf('errors:\n (FD) %.2e;\n (GB) %.2e;\n (LB) %.2e;\n (GP) %.2e;\n (WF) %.2e;\n',errz{:})
figure(7);clf
plot(x,Yp_fd,'r','linewidth',1);hold on
plot(x,Yp_glob,'y','linewidth',3)
plot(x,Ylp_loc,'g','linewidth',3)
plot(x,Yp_gp,'linewidth',3)
plot(x,Yp_wf,'c-.','linewidth',3)
plot(x,Yp,'k--','linewidth',3)
hold off
legend({'FD','GB','LP','GP','WF','true'})

%% functions

function [Globalderiv,Proj_mat,Phi_mat] = get_glob(p,x,varargin)
    inp = inputParser;
    addRequired(inp,'p');
    addRequired(inp,'x');
    addParameter(inp,'basis','poly');
    parse(inp,p,x,varargin{:})

    if isequal(inp.Results.basis,'poly')
        Phi = arrayfun(@(j)@(x)x.^j,0:p,'uni',0);
        Phip = arrayfun(@(j)@(x)j*x.^max(j-1,0),0:p,'uni',0);
    elseif isequal(inp.Results.basis,'sine')
        Phi = arrayfun(@(j)@(x)sin(j*x),1:p+1,'uni',0);
        Phip = arrayfun(@(j)@(x)j*cos(j*x),1:p+1,'uni',0);
    elseif isequal(inp.Results.basis,'cos')
        Phi = arrayfun(@(j)@(x)cos(j*x),0:p,'uni',0);
        Phip = arrayfun(@(j)@(x)-j*sin(j*x),0:p,'uni',0);
    elseif isequal(inp.Results.basis,'trig')
        Phi = [arrayfun(@(j)@(x)cos(j*x),0:floor(p/2),'uni',0),arrayfun(@(j)@(x)sin(j*x),1:ceil(p/2),'uni',0)];
        Phip = [arrayfun(@(j)@(x)-j*sin(j*x),0:floor(p/2),'uni',0),arrayfun(@(j)@(x)j*cos(j*x),1:ceil(p/2),'uni',0)];
    end
        
    %%% get proj mat
    Phi_mat = cell2mat(cellfun(@(phi)phi(x(:)),Phi,'uni',0));
    Phip_mat = cell2mat(cellfun(@(phi)phi(x(:)),Phip,'uni',0));
    Pinv = pinv(Phi_mat);
    Globalderiv = Phip_mat*Pinv;
    Proj_mat = Phi_mat*Pinv;

end

function Localderiv = get_loc(p,x,varargin)
    inp = inputParser;
    addRequired(inp,'p');
    addRequired(inp,'x');
    addParameter(inp,'basis','poly');
    parse(inp,p,x,varargin{:})

    if isequal(inp.Results.basis,'poly')
        Phi = arrayfun(@(j)@(x)x.^j,0:p,'uni',0);
        Phip = arrayfun(@(j)@(x)j*x.^max(j-1,0),0:p,'uni',0);
    elseif isequal(inp.Results.basis,'sine')
        Phi = arrayfun(@(j)@(x)sin(j*x),1:p+1,'uni',0);
        Phip = arrayfun(@(j)@(x)j*cos(j*x),1:p+1,'uni',0);
    elseif isequal(inp.Results.basis,'cos')
        Phi = arrayfun(@(j)@(x)cos(j*x),0:p,'uni',0);
        Phip = arrayfun(@(j)@(x)-j*sin(j*x),0:p,'uni',0);
    elseif isequal(inp.Results.basis,'trig')
        Phi = [arrayfun(@(j)@(x)cos(j*x),0:floor(p/2),'uni',0),arrayfun(@(j)@(x)sin(j*x),1:ceil(p/2),'uni',0)];
        Phip = [arrayfun(@(j)@(x)-j*sin(j*x),0:floor(p/2),'uni',0),arrayfun(@(j)@(x)j*cos(j*x),1:ceil(p/2),'uni',0)];
    end
    
    M = length(x);
    xlocal = x(1:p+1);
    xc = x(p/2+1);

    V = cell2mat(cellfun(@(phi)phi(xlocal(:)'),Phi(:),'uni',0));
    d = cellfun(@(phi)phi(xc),Phip(:));    
    W = V \ d;
    Localderiv = antisymconvmtx(flipud(W),M);
end

function [Mat,Proj_mat] = get_loc_lsq(p,x,M,varargin)
    m = floor((length(x)-1)/2);
    [Globalderiv,Proj_mat] = get_glob(p,x,varargin{:});
    Mat = antisymconvmtx(flipud(Globalderiv(m+1,:)'),M);
    Proj_mat = antisymconvmtx(flipud(Proj_mat(m+1,:)'),M);
end

function V = antisymmat(V,m)
    k = V(1,1:2*m+1);
    L = 2*m+1;
    V(:,m+2:L) = V(:,m+2:L) - fliplr(V(:,1:m));
    V(:,end-L+1:end-m-1) = V(:,end-L+1:end-m-1) - fliplr(V(:,end-m+1:end));
    V = V(:,m+1:end-m);
    lsum = flipud(2*cumsum(k(1:m))');
    rsum = 2*cumsum(fliplr(k(end-m+1:end)))';
    V(1:m,1) = V(1:m,1) + lsum;
    V(end-m+1:end,end) = V(end-m+1:end,end) + rsum;
end

function V = antisymconvmtx(k,N)
    V = sparseconvmat(k,N);
    L = length(k);
    m = (L-1)/2;
    V = antisymmat(V,m);
end

function V = sparseconvmat(k,N)
    k = flipud(k(:))';
    V = spdiags(repmat(k,N,1),0:length(k)-1,N,N+length(k)-1);
end

function V = symconvmtx(k,N)
    k = flipud(k(:))';
    L = length(k);
    l = (L-1)/2;
    V = spdiags(repmat(k,N,1),0:L-1,N,N+L-1);
    V(:,l+2:L) = V(:,l+2:L) + fliplr(V(:,1:l));
    V(:,end-L+1:end-l-1) = V(:,end-L+1:end-l-1) + fliplr(V(:,end-l+1:end));
    V = V(:,l+1:end-l);
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