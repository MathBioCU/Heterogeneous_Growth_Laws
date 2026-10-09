%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger


function tf = get_VVp_tf(Uobj,phifun,subinds_svd,rad_type,mt_params,lhs,K_min,toggle_VVp_svd)
    if isempty(lhs)
        lhs = eye(Uobj.nstates);
        lhs = mat2cell(lhs,ones(1,Uobj.nstates),Uobj.nstates);
    end

    numeq = length(lhs);

    tf = cellfun(@(ls)arrayfun(@(tf_param)testfcn(Uobj,'phifuns',phifun,'subinds',subinds_svd,...
        'meth',rad_type,'param',tf_param,'stateind',find(ls)),mt_params),lhs,'uni',0);
    mt = cell2mat(cellfun(@(t)[t.rads],tf(:),'uni',0));
    
    [cm,cn] = size(mt);
    V_cell = cell(numeq,1);
    Vp_cell = cell(numeq,1);
    for nn=1:numeq
        V_cell{nn} = [];
        Vp_cell{nn} = [];
        for j=1:length(mt_params)
            V_cell{nn} = [V_cell{nn};tf{nn}(j).get_testmat(0)];
        end
        [V_cell{nn},Vp_cell{nn}] = VVp_svd(full(V_cell{nn}),K_min,Uobj.grid{1},toggle_VVp_svd);
    end
    tf = cellfun(@(V,Vp)testfcn_Vcell(Uobj,{{V,-Vp}},{[0 1]}),V_cell,Vp_cell,'uni',0); 
end

function [V,Vp] = VVp_svd(V,K_min,t,toggle_VVp_svd)
    m = length(t);
    dt = mean(diff(t));
    [U,S,~] = svd(V',0);
    sings = diag(S);
    if toggle_VVp_svd>0
        s = find(cumsum(sings.^2)/sum(sings.^2)>toggle_VVp_svd^2,1);
        if isempty(s)
            s = min(K,size(V,1));
        end
    else
        corner_data = cumsum(sings)/sum(sings);
        s = getcorner_svd(corner_data,(1:length(corner_data))');%
        s = min(max(K_min,s),size(V,1));
    end
    inds = 1:s;
    V = U(:,inds)'*dt;

    Vp = V';
    Vp_hat = fft(Vp);
    if mod(m,2)==0
        k = [0:m/2 -m/2+1:-1]';
    else
        k = [0:floor(m/2) -floor(m/2):-1]';
    end
    Vp_hat = ((2*pi/m/dt)*1i*k).*Vp_hat;
    if mod(m,2)==0
        Vp_hat(m/2) = 0;        
    end
    Vp = real(ifft(Vp_hat))';

end

function tstarind = getcorner_svd(Ufft,xx)
    NN = length(Ufft);
    Ufft = Ufft/max(abs(Ufft))*NN;
    errs = zeros(NN,1);
    for k=1:NN
        [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k);
        errs(k) = sum(abs((L1-Ufft_av1)./Ufft_av1)) + sum(abs((L2-Ufft_av2)./Ufft_av2)); % relative l1           
    end
    [~,tstarind] = min(errs);
end

function [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k)
   NN = length(Ufft);
   subinds1 = 1:k;
   subinds2 = k:NN;
   Ufft_av1 = Ufft(subinds1);
   Ufft_av2 = Ufft(subinds2);
   
   [m1,b1,L1]=lin_regress(Ufft_av1,xx(subinds1));
   [m2,b2,L2]=lin_regress(Ufft_av2,xx(subinds2));
end

function [m,b,L]=lin_regress(U,x)
   m = (U(end)-U(1))/(x(end)-x(1));
   b = U(1)-m*x(1);
   L = U(1)+m*(x-x(1));
end