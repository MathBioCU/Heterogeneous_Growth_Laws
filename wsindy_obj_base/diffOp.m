classdef diffOp < linearOp
    properties
        difftags
        meth
        stateind
        Dmats
        w
        ndims
    end

    methods
        function obj = diffOp(difftags,varargin)
            default_nstates = 1;
            default_stateind = 1;
            default_meth = 'fd';
            default_gradterms = {};
            default_gradon = 0;
            default_w = [];
            default_dat = [];

            p = inputParser;
            addRequired(p,'difftags');
            addParameter(p,'nstates',default_nstates);
            addParameter(p,'stateind',default_stateind);
            addParameter(p,'meth',default_meth);
            addParameter(p,'gradterms',default_gradterms);
            addParameter(p,'gradon',default_gradon);
            addParameter(p,'w',default_w);
            addParameter(p,'dat',default_dat);
            parse(p,difftags,varargin{:})

            obj.Dmats = [];
            obj.difftags = difftags;
            obj.ndims = length(difftags);
            obj.stateind = p.Results.stateind;
            obj.nstates = max(obj.stateind,p.Results.nstates);
            obj.fHandle = @(varargin)varargin{obj.stateind};
            obj.gradterms = p.Results.gradterms;
            obj.meth = p.Results.meth;
            obj.gradon = p.Results.gradon;
            obj.w = p.Results.w;
            dat = p.Results.dat;
            if ~isempty(dat)
                obj = obj.get_Dmats(dat);
            end

            if and(obj.gradon,isempty(obj.gradterms))
                obj = obj.get_grads;
            end
            obj.linOp = obj;
        end

    end

    methods

        function Y = evalterm(obj,dat,varargin)

            if isequal(class(dat),'cell')
                Y = dat{obj.stateind};
            elseif isequal(class(dat),'wsindy_data')
                Y = dat.Uobs{obj.stateind};            
            else
                Y = dat(:,obj.stateind);
            end

            if isempty(obj.Dmats)
                default_get_Dmats = 1;
            elseif size(obj.Dmats{1},2)~=size(Y,1)
                default_get_Dmats = 1;
            else
                default_get_Dmats = 0;
            end

            p = inputParser;
            addRequired(p,'dat');
            addParameter(p,'get_Dmats', default_get_Dmats);
            parse(p,dat,varargin{:})

            if p.Results.get_Dmats
                obj = obj.get_Dmats(dat);
            end

            for j=1:length(obj.Dmats)
                if obj.difftags(j)~=0
                    shift = 1:length(obj.Dmats);
                    shift([1 j]) = [j 1];
                    Y = permute(Y,shift);
                    Y = pagemtimes(full(obj.Dmats{j}),Y);
                    Y = permute(Y,shift);
                end
            end
        end

        function s = get_str(obj)
            s = ['(d/dt)^',num2str(obj.difftags),'(x',num2str(obj.stateind),')'];
            % 
            % if obj.stateind==1
            %     s = ['(d/dt)^',num2str(obj.difftags),'(x)'];
            % elseif obj.stateind==2
            %     s = ['(d/dt)^',num2str(obj.difftags),'(z)'];
            % end

        end

        function Y = evaltermLinOp(obj,dat)
            Y = obj.evalterm(dat);
        end

        function Y = evalgrads(obj,dat)
            if isempty(obj.gradterms)
                obj.get_grads;
                obj.gradon = 1;
            end
            if isequal(class(dat),'wsindy_data')
                e = zeros(1,nstates);
                e(obj.stateind) = 1;
                e = num2cell(e);
                Y = cellfun(@(u,e)u*e,dat.Uobs,e,'uni',0);
            else 
                Y = dat*0;
                n = find(size(Y)==obj.nstates);
                if n(1)==1
                    Y(obj.stateind,:)=1;
                    Y = mat2cell(Y,ones(obj.nstates,1),size(Y,2))';
                elseif n(1)==2
                    Y(:,obj.stateind)=1;
                    Y = mat2cell(Y,size(Y,1),ones(1,obj.nstates));
                end
            end
            % dims = num2cell(dat.dims);
            % Y = repmat({sparse(dims{:})},1,dat.nstates);
            % Y{obj.stateind} = obj.evalterm(dat);
        end

        function obj = get_Dmats(obj,dat)
            if isnan(obj.difftags)
                obj.Dmats = arrayfun(@(i)sparse(dat.dims(i),dat.dims(i)),(1:dat.ndims)','uni',0);
            else   
                if isequal(obj.meth,'fd')
                    obj.Dmats = obj.fdMats(dat);
                elseif isequal(obj.meth,'wf')
                    obj.Dmats = obj.wfMats(dat);
                elseif isequal(obj.meth,'wffd')
                    obj.Dmats = obj.wffdMats(dat);
                end
            end
        end

        function Y = diffmat(obj,dat)
            if isempty(obj.Dmats)
                obj = obj.get_Dmats(dat);
            end
            Y = repmat({sparse(prod(dat.dims),prod(dat.dims))},1,obj.nstates);
            Z = 1;
            for k=1:dat.ndims
                Z = kron(obj.Dmats{k},Z);
            end
            Y{obj.stateind} = Z;
        end

        function rhs = get_rhs(obj)
            if length(obj.difftags)==1
                x = obj.stateind + obj.difftags*obj.nstates;
                rhs = @(varargin) varargin{x};
            end
        end


    end

    methods (Hidden = true)

        function [Dmats,c] = fdMats(obj,dat,wd)
            if ~exist('wd','var')
                if isempty(obj.w)
                    obj.w = max(obj.difftags);
                end
                wd = obj.w;
            end
            if all(dat.isUniform)
                Dmats = cell(dat.ndims,1);
                c = cell(dat.ndims,1);
                for i=1:dat.ndims
                    dv = mean(diff(dat.grid{i}));
                    c{i} = fdcoeffF(obj.difftags(i),0,(-wd:wd)*dv);
                    Dmats{i} = obj.antisymconvmtx(flipud(c{i}(:,end)),dat.dims(i));
                end
            end
        end

        function Dmats = wfMats(obj,dat)
            if isequal(class(dat),'double')
                sz = size(dat);
                ndims = length(sz);
                grids = {};
                for j=1:ndims
                    if sz(j)>1
                        grids = [grids,{(0:sz(j)-1)'}];
                    end
                end
                dat = wsindy_data(dat,grids);
            end
            if all(dat.isUniform)
                Dmats = cell(dat.ndims,1);
                dat.getcorners;
                for i=1:dat.ndims
                    k_x = dat.ks(i);
                    dv = mean(diff(dat.grid{i}));

                    eta = 9;
                    phi = @(t) (1-t.^2).^eta;
                    phip = @(t) (-2*eta)*(t.*(1-t.^2).^(eta-1));
                    if isempty(obj.w)
                        tauhat = 2;
                        % obj.w = get_tf_support(phi, dat.dims(i), tauhat, k_x);
                        obj.w = ceil(tauhat*dat.dims(i)*sqrt(eta*2+3)/2/pi/k_x);
                    end
                    xf = linspace(-1,1,2*obj.w+1);
                    Cfs = zeros(2,2*obj.w+1);
                    Cfs(1,:) = phi(xf);
                    if obj.difftags(i)>1
                        syms y;
                        phip = matlabFunction(diff(phi(y),obj.difftags(i)));
                    end
                    Cfs(2,:) = phip(xf)*(obj.w*dv).^-obj.difftags(i);

                    Vp = obj.antisymconvmtx(Cfs(2,:),dat.dims(i));
                    V = obj.symconvmtx(Cfs(1,:),dat.dims(i));

                    [Uu,Ss,Vv] = svd(full(V),0);
                    % figure
                    % findchangepts(diag(Ss),'Statistic','linear')
                    ss = findchangepts(diag(Ss),'Statistic','linear');
                    % disp(norm(diag(Ss(1:ss,1:ss)))/norm(diag(Ss)))
                    % ss = size(Vv,2);
                    Dmats{i} = Vv(:,1:ss)*diag(1./diag(Ss(1:ss,1:ss)))*Uu(:,1:ss)'*Vp;

                                        % s = ceil(obj.w/5);
                    % Vp = Vp(1:s:end,:); V = V(1:s:end,:);
                    % plot(diag(Ss))
                    % ss = getcorner(diag(Ss),1:size(Ss,2));
                    % F = dftmtx(dat.dims(i));
                    % A = real(F(:,1:k_x)*F(:,1:k_x)'*Vv*Ss)/dat.dims(i);
                    % size(A)
                    % plot(fliplr(vecnorm(A-Vv*Ss)./vecnorm(Vv*Ss)))
                    % ss = getcorner(fliplr(vecnorm(A-Vv*Ss)),1:size(A,2));
                    % ss = size(A,2) - ss/2 + 1;
                    % ss = floor(ss);
                    bc_ind = 3;%ceil(obj.w/2);
                    bc = [repmat(Dmats{i}(bc_ind+1,:),bc_ind,1);...
                        repmat(Dmats{i}(dat.dims(i)-bc_ind,:),bc_ind,1)];
                    Dmats{i}([1:bc_ind dat.dims(i)-bc_ind+1:dat.dims(i)],:) = bc;

                end
            end
        end

        function Dmats = wffdMats(obj,dat)
            if isequal(class(dat),'double')
                sz = size(dat);
                ndims = length(sz);
                grids = {};
                for j=1:ndims
                    if sz(j)>1
                        grids = [grids,{(0:sz(j)-1)'}];
                    end
                end
                dat = wsindy_data(dat,grids);
            end
            if all(dat.isUniform)

                p_loc = 6; % must be even
                [Dmats_fd,c_fd] = obj.fdMats(dat,p_loc);
                Dmats = cell(dat.ndims,1);
                dat.getcorners;
                sigest = dat.estimate_sigma;
                for i=1:dat.ndims
                    if obj.difftags~=0
                        k_x = dat.ks(i);
                        dv = mean(diff(dat.grid{i}));
                        [phi,phip] = optTFcos(3,0);
                        if isempty(obj.w)
                            tauhat = 1+dat.dims(i)^(1/3)*sigest{obj.stateind};
                            obj.w = get_tf_support(phi, dat.dims(i), tauhat, k_x);
                            obj.w = max(4,obj.w);
                        end
                        xf = linspace(-1,1,2*obj.w+1);
                        Cfs = zeros(2,2*obj.w+1);
                        Cfs(1,:) = phi(xf);
                        if obj.difftags(i)>1
                            syms y;
                            phip = matlabFunction(diff(phi(y),obj.difftags(i)));
                        end
                        Cfs(2,:) = phip(xf)*(obj.w*dv).^-obj.difftags(i);
                        Vp = obj.antisymconvmtx(Cfs(2,:)/norm(Cfs(1,:),1),dat.dims(i));
                        if mod(obj.difftags(i),2)==1
                            V = obj.symconvmtx(Cfs(1,:)/norm(Cfs(1,:),1),dat.dims(i));
                        else
                            V = obj.antisymconvmtx(Cfs(1,:)/norm(Cfs(1,:),1),dat.dims(i));
                        end
                        c = 0.00001/(0.1+sigest{obj.stateind});
                        Localderiv = obj.antisymconvmtx(flipud(c_fd{i}(:,2)),dat.dims(i));
                        regmat1 = (max(sigest{obj.stateind},10^-6)/norm(Localderiv(1,:),1))*Localderiv;
                        regmat2 = c*speye(dat.dims(i));
                        Bmat1 = sparse(dat.dims(i),dat.dims(i));
                        Bmat2 = c*Dmats_fd{i};
                        A = [V;regmat1;regmat2];
                        B = [Vp;Bmat1;Bmat2];
                        Dmats{i} = A \ B;
                    else
                        Dmats{i} = Dmats_fd{i};
                        V = speye(dat.dims(i));             
                        Vp = speye(dat.dims(i));
                    end
                    c = 0.00001/(0.1+sigest{obj.stateind});
                    Localderiv = obj.antisymconvmtx(flipud(c_fd{i}(:,end)),dat.dims(i));
                    regmat1 = (max(sigest{obj.stateind},10^-6)/norm(Localderiv(1,:),1))*Localderiv;
                    regmat2 = c*speye(dat.dims(i));
                    Bmat1 = sparse(dat.dims(i),dat.dims(i));
                    Bmat2 = c*Dmats_fd{i};
                    A = [V;regmat1;regmat2];
                    B = [Vp;Bmat1;Bmat2];
                    Dmats{i} = A \ B;
                end
            end
        end

        function V = symconvmtx(obj,k,N)
            k = flipud(k(:))';
            L = length(k);
            l = (L-1)/2;
            V = spdiags(repmat(k,N,1),0:2*obj.w,N,N+2*obj.w);
            V(:,l+2:L) = V(:,l+2:L) + fliplr(V(:,1:l));
            V(:,end-L+1:end-l-1) = V(:,end-L+1:end-l-1) + fliplr(V(:,end-l+1:end));
            V = V(:,l+1:end-l);
        
        end

        function V = antisymconvmtx(obj,k,N)
            k = flipud(k(:))';
            L = length(k);
            l = (L-1)/2;
            V = spdiags(repmat(k,N,1),0:L-1,N,N+L-1);
            V(:,l+2:L) = V(:,l+2:L) - fliplr(V(:,1:l));
            V(:,end-L+1:end-l-1) = V(:,end-L+1:end-l-1) - fliplr(V(:,end-l+1:end));
            V = V(:,l+1:end-l);
            lsum = flipud(2*cumsum(k(1:l))');
            rsum = (2*cumsum(fliplr(k(end-l+1:end)))');
            V(1:l,1) = V(1:l,1) + lsum;
            V(end-l+1:end,end) = V(end-l+1:end,end) + rsum;
        end

        function obj = get_grads(obj)
            obj.gradterms = repelem(diffOp(NaN,'nstates',obj.nstates),1,obj.nstates);
            obj.gradterms(obj.stateind) = obj;
            obj.gradon = 1;
        end

        function Cfs = phi_weights(obj,k,diffs)
            xf = linspace(-1,1,2*obj.rads(k)+1);
            x = xf(2:end-1);
            Cfs = zeros(length(diffs),2*obj.rads(k)+1);
            syms y;
            f = @(y)obj.phifuns{k}(y);
            for j=1:length(diffs)
                Df = matlabFunction(diff(f(y),diffs(j)));
                Cfs(j,2:end-1) = fillmissing(Df(x),'constant',Df(eps));
                inds = find(isinf(abs(Cfs(j,:))));
                for i=1:length(inds)
                    Cfs(j,inds(i)) = Df(xf(inds(i))-sign(xf(inds(i)))*eps);
                end
            end
            Cfs = Cfs.*(obj.rads(k)*obj.dv(k)).^(-diffs(:));
%             Cfs = Cfs/norm(Cfs(1,:),1);
        end

        function Localderiv = get_loc(obj,p,x,varargin)
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
        
            %%% get proj mat
            V = cell2mat(cellfun(@(phi)phi(xlocal(:)'),Phi(:),'uni',0));
            d = cellfun(@(phi)phi(xc),Phip(:));    
            W = V \ d;
            Localderiv = obj.antisymconvmtx(flipud(W),M);
            % Localderiv = spdiags(repmat(w(:)',M,1),[0:p],M,M+p);
        end

        function m = get_scale(obj,scales)
            m = 1;
            if ~isempty(scales)
                m=m*scales(obj.stateind);
                m=m/prod(scales(obj.nstates+1:end).^obj.difftags);
            end
        end

    end

end
