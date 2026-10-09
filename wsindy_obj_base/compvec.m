classdef compvec < term
    properties
        t1
        t2
    end

    methods
        function obj = compvec(t2,t1,varargin)
            % written as t2(t1(x))
            obj = obj@term('gradon',0);
            default_gradon = 0;
            default_gradterms = {};
            p = inputParser;
            addRequired(p,'t1');
            addRequired(p,'t2');
            addParameter(p,'gradterms',default_gradterms);
            addParameter(p,'gradon',default_gradon);
            parse(p,t1,t2,varargin{:})

            if isequal(class(t2),'function_handle')
                t2 = term('fHandle',t2);
            elseif isequal(class(t2),'double')
                t2 = term('ftag',t2);
            end

            if length(t1.terms)~=t2.nstates
                disp('unable to match functions')
                return
            end

            %%% class(t1) must be 'library'
            obj.t1 = t1;
            obj.t2 = t2;
            obj.nstates = t1.nstates;
            obj.gradterms = p.Results.gradterms;
            obj.gradon = p.Results.gradon;
            obj.linOp = t2.linOp;
            obj.t2.linOp = [];
            obj.ftag = {'compvec',t2.ftag};
            obj.fHandle = @(dat)obj.evalterm(dat);

            if and(obj.gradon,isempty(obj.gradterms))
                obj = obj.get_grads;
            end
        end
    end

    methods

        function Y = evalterm(obj,dat)
            Y = obj.t1.evalterms(dat);
            Y = mat2cell(Y,size(Y,1),ones(1,size(Y,2)));
            Y = obj.t2.evalterm(Y);
        end

        function rhs = get_rhs(obj)
            rhs1 = cellfun(@(tt)tt.get_rhs,obj.t1,'uni',0);
            rhs = @(varargin)obj.get_rhs_in(rhs1,obj.t2.fHandle,varargin{:});
            % rhs1 = @(varargin)cellfun(@(tt)tt(varargin{:}),rhs1);
            % rhs = @(varargin)obj.t2.fHandle(rhs1(varargin{:}));
        end

        function out = get_rhs_in(obj,fcell,flast,varargin)
            out = cellfun(@(f)f(varargin{:}),fcell,'uni',0);
            out = flast(out{:});
        end


    end

end
