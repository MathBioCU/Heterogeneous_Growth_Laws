classdef compterm < term
    properties
        t1
        t2
    end

    methods
        function obj = compterm(t2,t1,varargin)
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

            if isequal(class(t1),'function_handle')
                t1 = term('fHandle',t1);
            elseif isequal(class(t1),'double')
                t1 = term('ftag',t1);
            end
            if isequal(class(t2),'function_handle')
                t2 = term('fHandle',t2);
            elseif isequal(class(t2),'double')
                t2 = term('ftag',t2);
            end
            obj.t1 = t1;
            obj.t2 = t2;
            obj.nstates = t1.nstates;
            obj.gradterms = p.Results.gradterms;
            obj.gradon = p.Results.gradon;
            obj.linOp = t2.linOp;
            obj.t2.linOp = [];
            obj.ftag = {'comp',t1.ftag,t2.ftag};
            obj.fHandle = @(dat)obj.evalterm(dat);

            if and(obj.gradon,isempty(obj.gradterms))
                obj = obj.get_grads;
            end
        end
    end

    methods

        function Y = evalterm(obj,dat)
            if isequal(class(obj.t2),'diffOp')
                Y = obj.t1.evaltermLinOp(dat);
            else
                Y = obj.t2.evalterm(obj.t1.evaltermLinOp(dat));
            end
        end

        function obj = get_grads(obj)
            obj.gradterms = repmat(prodterm.empty,1,obj.nstates);
            if obj.t1.gradon==0
                obj.t1.get_grads;
            end
            if obj.t2.gradon==0
                obj.t2.get_grads;
            end
            for j=1:obj.nstates
                obj.gradterms(j) = prodterm(compterm(obj.t2.gradterms,obj.t1),obj.t1.gradterms(j),'gradon',0);
            end
        end

        function s = get_str(obj)
            s = ['(',obj.t2.get_str,')o(',obj.t1.get_str,')'];
            if isequal(class(obj.linOp),'diffOp')
                s = ['(d/dt)^',num2str(obj.linOp.difftags),s];
            end
        end

        function Y = diffmat(obj,dat)
            if isequal(class(obj.t2),'diffOp')
                Y=obj.t1.diffmat(dat);
            else
                Z=obj.t2.diffmat(obj.t1.evalterm(dat));
                Y=obj.t1.diffmat(dat);
                Y=cellfun(@(L) Z{1}*L,Y,'uni',0);
            end
        end

        function rhs = get_rhs(obj)
            rhs1 = obj.t1.get_rhs;
            rhs = @(varargin)obj.t2.fHandle(rhs1(varargin{:}));
        end

        function m = get_scale(obj,scales)
            m1 = obj.t1.get_scale(scales);
            m2 = obj.t2.get_scale(scales);
            if m2 ~= 1
                m = m2*m1;
            else
                m = 1;
            end
        end
        
    end

end
