% tags = get_tags(1:2,[],4)

function tags = get_tags(polys,trigs,n,varargin)
    p = inputParser;
    addRequired(p,'polys');
    addRequired(p,'trigs');
    addRequired(p,'n');
    addParameter(p,'uni',1);
    addParameter(p,'neg',0);
    addParameter(p,'boolT',@(tags)repmat(1>0,size(tags,1),1));
    parse(p,polys,trigs,n,varargin{:})
    uni = p.Results.uni;
    neg = p.Results.neg;
    boolT= p.Results.boolT;

    tags = [];
    for p = 1:length(polys)
        monom_powers = partitionNk(polys(p),n);
        tags = [tags;monom_powers];
    end

    if neg==1
        a = (-ones(1,n)).^de2bi(0:2^n-1);
        tags = cell2mat(arrayfun(@(i)tags*diag(a(i,:)),(1:2^n)','uni',0));
    end
    tags = unique(tags(boolT(tags),:),'rows');

    if isvector(trigs)
        for k=1:length(trigs)
            trig_inds = [-trigs(k)*1i*eye(n);trigs(k)*1i*eye(n)];
            tags = [tags; trig_inds];
        end
        tags = unique(tags,'rows');
    else
        for k=1:size(trigs,1)
            trig_inds = [-1i*diag(trigs(k,:));1i*diag(trigs(k,:))];
            tags = [tags; trig_inds];
        end
        tags = unique(tags,'rows');
    end

    if ~uni
        tags = mat2cell(tags,ones(size(tags,1),1),size(tags,2));
    end

end
