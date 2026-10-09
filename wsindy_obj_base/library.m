%%% defines a library object simply used to gather a set of terms. Operations
%%% consist mainly of mapping term methods over a list of terms, with some exceptions
%%% (such as trimlib) which performs a library trimming based on taylor expansion of the terms

classdef library < handle

    properties
        terms
        tags
        nstates
    end

    methods
        
        function obj = library(varargin)
            p = inputParser;
            addParameter(p,'tags',{});
            addParameter(p,'nstates',1);
            addParameter(p,'terms',{});
            addParameter(p,'polys',[]);
            addParameter(p,'gradon',1);
            addParameter(p,'neg',0);
            addParameter(p,'boolT',@(tags)true(size(tags,1),1));
            parse(p,varargin{:})
            obj.tags = {};
            obj.terms = {};
            if ~isempty(p.Results.terms)
                obj = add_terms(obj,p.Results.terms,p.Results.gradon);
            end

            if ~isempty(p.Results.tags)
                obj = add_terms(obj,p.Results.tags,p.Results.gradon);
            end

            if ~isempty(obj.terms)
                obj.nstates = obj.terms{1}.nstates;
            else
                obj.nstates = p.Results.nstates;
            end

            if ~isempty(p.Results.polys)
                tags_poly = get_tags(p.Results.polys,[],obj.nstates,'neg',p.Results.neg,'boolT',p.Results.boolT);
                obj = add_terms(obj,tags_poly,p.Results.gradon);
            end

        end

    end

    methods
        function fs = get_fHandles(obj,inds)
            fs = cellfun(@(t)t.fHandle,obj.terms(inds),'uni',0);
        end

        function obj = add_tags(obj,varargin)

            p = inputParser;
            addParameter(p,'polys',[]);
            addParameter(p,'trigs',[]);
            addParameter(p,'uni',1);
            addParameter(p,'neg',0);
            addParameter(p,'boolT',[]);
            addParameter(p,'boolTL',[]);
            addParameter(p,'lhs',[]);
            addParameter(p,'gradOn',0)
            parse(p,varargin{:})
            polys = p.Results.polys;
            trigs = p.Results.trigs;
            uni = p.Results.uni;
            neg = p.Results.neg;
            boolT= p.Results.boolT;
            boolTL = p.Results.boolTL;
            lhs = p.Results.lhs;
            gradOn = p.Results.gradOn;
        
            tags_in = [];
            for p = 1:length(polys)
                monom_powers = partitionNk(polys(p),obj.nstates);
                tags_in = [tags_in;monom_powers];
            end
        
            if neg==1
                a = (-ones(1,obj.nstates)).^de2bi(0:2^obj.nstates-1);
                tags_in = cell2mat(arrayfun(@(i)tags_in*diag(a(i,:)),(1:2^obj.nstates)','uni',0));
            end

            if isvector(trigs)
                for k=1:length(trigs)
                    trig_inds = [-trigs(k)*1i*eye(obj.nstates);trigs(k)*1i*eye(obj.nstates)];
                    tags_in = [tags_in; trig_inds];
                end
                tags_in = unique(tags_in,'rows');
            else
                for k=1:size(trigs,1)
                    trig_inds = [-1i*diag(trigs(k,:));1i*diag(trigs(k,:))];
                    tags_in = [tags_in; trig_inds];
                end
                tags_in = unique(tags_in,'rows');
            end

            if ~isempty(boolT)
                tags_in = unique(tags_in(boolT(tags_in),:),'rows');
            end
            if ~isempty(boolTL)
                tags_in = unique(tags_in(boolTL(tags_in,lhs),:),'rows');
            end
        
            if ~uni
                tags_in = mat2cell(tags_in,ones(size(tags_in,1),1),size(tags_in,2));
            end
            obj.add_terms(tags_in,gradOn);

        end

        function obj = add_terms(obj,terms_in,gradon,dup)
            if ~exist('gradon','var')
                gradon = 1;
            end
            if ~exist('dup','var')
                dup = 0;
            end
            if isequal(class(terms_in),'double') % if collection of pw terms, then add
                terms_in = mat2cell(terms_in,ones(size(terms_in,1),1),size(terms_in,2));
            elseif isequal(class(terms_in),'function_handle')
                terms_in = term('fHandle',terms_in);
            end
            if isequal(class(terms_in),'cell')
                for j=1:length(terms_in)
                    if or(dup,~any(cellfun(@(t)isequal(t,terms_in{j}),obj.tags)))
                        if isequal(class(terms_in{j}),'double')
                            if length(terms_in{j})>obj.nstates
                                obj.terms = [obj.terms,{term('ftag',terms_in{j}(1:obj.nstates),'linOp',terms_in{j}(obj.nstates+1:end),'gradon',gradon)}];
                            else
                                obj.terms = [obj.terms,{term('ftag',terms_in{j},'gradon',gradon)}];
                            end
                            
                            obj.tags = [obj.tags,terms_in(j)];
                    elseif isequal(class(terms_in{j}),'function_handle')
                            obj.terms = [obj.terms,{term('fHandle',terms_in{j},'gradon',gradon)}];
                            obj.tags = [obj.tags,terms_in(j)];
                        elseif any(cellfun(@(x) isequal(x,'absterm'),superclasses(terms_in{j})))
                            obj.terms = [obj.terms,{terms_in{j}}];
                            obj.tags = [obj.tags,terms_in(j)];
                        end
                    end
                end
            elseif any(cellfun(@(x) isequal(x,'absterm'),superclasses(terms_in)))
                for j=1:length(terms_in)
                    obj.terms = [obj.terms,{terms_in(j)}];
                    if isequal(class(terms_in(j).linOp),'diffOp')
                        tag_in = [terms_in(j).ftag terms_in(j).linOp.difftags];
                    else
                        tag_in = terms_in(j).ftag;
                    end
                    obj.tags = [obj.tags,{tag_in}];
                end
            end
            if ~isempty(terms_in)
                obj.nstates = obj.terms{1}.nstates;
            end
        end

        function Theta = evalterms(obj,dat,S)
            if ~exist('S','var')
                S = true(length(obj.terms),1);
            end
            Theta = cell2mat(cellfun(@(tm) tm.evalterm(dat),obj.terms(S),'uni',0));
        end

        function obj = remove_terms(obj,ndims,varargin)
            p = inputParser;
            addParameter(p,'ftag',[]);
            addParameter(p,'difftag',[]);
            addParameter(p,'ftag_difftag',{});
            parse(p,varargin{:})



            ftag = p.Results.ftag;
            difftag = p.Results.difftag;
            ftag_difftag = p.Results.ftag_difftag;

            remove_inds = false(length(obj.terms),1);
            for k=1:length(obj.terms)
                dt = zeros(1,ndims);
                ft = obj.terms{k}.ftag;
                if isequal(class(obj.terms{k}.linOp),'diffOp')
                    dt = obj.terms{k}.linOp.difftags;
                end
                if ~isempty(ftag)
                    if any(cellfun(@(c)c(ft),ftag))
                        remove_inds(k) = true;
                    end
                end
                if ~isempty(difftag)
                    if any(cellfun(@(c)c(dt),difftag))
                        remove_inds(k) = true;
                    end
                end                
                if ~isempty(ftag_difftag)
                    if isequal(dt,ftag_difftag{1})
                        if any(cellfun(@(c)c(ft),ftag_difftag{2}))
                            remove_inds(k) = true;
                        end
                    end
                end                    
            end
            obj.terms = obj.terms(~remove_inds);
            obj.tags = obj.tags(~remove_inds);            
        end


        function Theta_grad = evalGradterms(obj,dat,S)
            if ~exist('S','var')
                S = true(length(obj.terms),1);
            end
            Theta_grad = cellfun(@(tm) tm.evalgrads(dat),obj.terms(S),'uni',0);
        end

        function Theta_grad = grad2mat(obj,dat,S)
            if ~exist('S','var')
                S = true(length(obj.terms),1);
            end
            Theta_grad = obj.evalGradterms(dat,S);
            Theta_grad = cellfun(@(tm) reshape([tm{:}],[],1), Theta_grad,'uni',0);
            Theta_grad = cell2mat(Theta_grad);
        end

        function Theta_H = evalHvec(obj,dat,J)
            Theta_grad = grad2mat(obj,dat);
            Theta_H = J*Theta_grad;
        end

        function obj = complib(obj,t2s,t1s)
            for j=1:length(t2s)
                for i=1:length(t1s)
                    obj.add_terms(compterm(t2s{j},t1s{i}));
                end
            end
        end

        function obj = trimlib(obj,pt,ord,tol,toggleH,togglepar)
            if ~exist('toggleH','var')
                toggleH=0;
            end
%             q = sym('q');
            args = num2cell(pt);
            tays = cell(1,length(obj.tags));
            ns = obj.nstates;
            termz = obj.terms;
            parfor(j=1:length(obj.tags),togglepar)
                tays{j} = termz{j}.get_taylor(pt,ord,toggleH);
            end
            tays = cell2mat(tays);
            sp = tays~=0;
            Xsparse = unique(sp','rows')';
            ii = cell(size(Xsparse,2),1);
            As = {};
            if toggleH==1
                ord = ord -1;
            end
            for j=1:size(Xsparse,2)
                ii = find(~any(sp-Xsparse(:,j),1));
                if length(ii)>1
                    G = tays(:,ii);
                    for l=1:obj.nstates
                        k = find(Xsparse((l-1)*ord+(1:ord),j),1);
                        if ~isempty(k)
                            G((l-1)*ord+(1:ord),:) = G((l-1)*ord+(1:ord),:)./G((l-1)*ord+k,:);
                        end
                    end
                    %%% compute ||G_i-G_j||, note the first non-zero entries already agree
                    A = real(sqrt((vecnorm(G)').^2 - 2*(G')*G + (vecnorm(G).^2)));
                    A = A-triu(A);
                    S = ones(size(A));
                    S = logical(S - triu(S));
                    A = and(A < tol,S);
                    if any(A(:))
                        As = [As,{{Xsparse(:,j),ii,G,A}}];
                    end
                end
            end
            
            rm_inds = [];
            for i=1:length(As)
                max_inds = [];
                for j=1:obj.nstates
                    foo_ind = 1+ord-find(flipud(As{i}{1}((j-1)*ord+1:j*ord)),1);
                    max_inds = [max_inds foo_ind];
                end
                while any(As{i}{end}(:))
                    r={[],[]};
                    [r{:}] = find(As{i}{end});
                    r = [r{:}];
                    r = unique(r(:));
                    [~,a] = max(sum(abs(As{i}{3}(max_inds,r))));
                    rm = r(a);
                    rm_inds = [rm_inds As{i}{2}(rm)];
                    As{i}{2} = As{i}{2}([1:rm-1 rm+1:end]);
                    As{i}{3} = As{i}{3}(:,[1:rm-1 rm+1:end]);
                    As{i}{4} = As{i}{4}(:,[1:rm-1 rm+1:end]);
                    As{i}{4} = As{i}{4}([1:rm-1 rm+1:end],:);
                end
            end
            
            J = length(obj.terms);
            obj.terms = obj.terms(setdiff(1:J,rm_inds)); 
            obj.tags = obj.tags(setdiff(1:J,rm_inds)); 
        end

        function Mscales = get_scales(obj,scales)
            Mscales = cellfun(@(t)t.get_scale(scales),obj.terms(:));
        end

        function diff_tags = get_diffs(obj,ndims)
            diff_tags = [];
            for tm = [obj.terms{:}]
                if isequal(class(tm.linOp),'diffOp')
                    diff_tags = [diff_tags;tm.linOp.difftags];
                elseif isempty(tm.linOp)
                    diff_tags = [diff_tags;zeros(1,ndims)];
                end
            end
            diff_tags = unique(diff_tags,'rows');
        end
    end

end