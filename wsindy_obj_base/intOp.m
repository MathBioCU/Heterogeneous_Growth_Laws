classdef intOp < linearOp
    properties
        f_in % inner function of state
        K % function of time
        Kgrid
        ndims
        Kmat
    end

    methods
        function obj = intOp(ndims,varargin)
            p = inputParser;
            addRequired(p,'ndims');
            addParameter(p,'K',[]);
            addParameter(p,'nstates',1);
            addParameter(p,'f_in',[]);
            parse(p,ndims,varargin{:})

            obj.nstates = p.Results.nstates;
            obj.ndims = p.Results.ndims;

            obj.f_in = p.Results.f_in;
            if isempty(obj.f_in)
                obj.f_in = term(ftag,zeros(1,obj.nstates));
            end

            obj.K = p.Results.K;
            if isempty(obj.K)
                obj.K = @(varargin) varargin{1}*0;
            end
            obj.Kmat = [];
            obj.Kgrid = [];
        end

    end

    methods

        function Y = evalterm(obj,dat)
            if dat.dims~=length(obj.Kmat)
                obj.get_Kmat(dat);
            end            
            Y = obj.f_in.evalterm(dat);
            Y = obj.Kmat*Y(:);
        end

        function get_Kmat(obj,dat)
            obj.Kgrid = obj.K(dat.grid{1});
            K_tri = zeros(dat.dims, dat.dims);
            for i=1:dat.dims
                K_tri(i,1:i) = obj.Kgrid(i:-1:1)*dat.dv;
                K_tri(i,[1 i]) = K_tri(i,[1 i])/2;
            end
            obj.Kmat = K_tri;
        end

        function s = get_str(obj)
            f_in_tag = functions(obj.f_in.fHandle).function;
            ind = strfind(f_in_tag,')');
            s = [strrep(functions(obj.K).function,'@(t)',''),'*',f_in_tag(ind(1)+1:end)];
        end

        function plot_kernel(obj,dat)
            plot(dat.grid{1},obj.Kgrid)
        end

        function m = get_scale(obj,scales)
            m = 1;
            if ~isempty(scales)
                m = obj.f_in.get_scale(scales);
                m = m*scales(end);
            end
        end

        function [rhs,t,Kweights] = get_rhs(obj,tau,N)
            t = (0:tau:(N-1)*tau)';
            Kweights = obj.K(t)*tau;
            Kweights([1,end]) = Kweights([1,end])/2;
            rhs = @(y) dot(obj.f_in.evalterm(y),Kweights);            
        end
    end

end
            % ndg  = cell(obj.ndims,1);
            % [ndg{:}] = ndgrid(dat.grid{:});
            % Kgrid = obj.K(ndg{:});
            
            % Y = convn(Y,Kgrid,'same');


