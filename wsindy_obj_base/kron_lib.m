function L = kron_lib(L1,L2)
    [a,b] = ndgrid(1:length(L1.terms),1:length(L2.terms));
    L = library();
    arrayfun(@(i,j)L.add_terms(kron_terms(L1.terms{i},L2.terms{j})),a(:),b(:),'uni',0);
end

function tout = kron_terms(t1,t2)
    if all([~isempty(t1.ftag) isnumeric(t1.ftag) ~isempty(t2.ftag) isnumeric(t2.ftag)])
        tout = term('ftag',[t1.ftag t2.ftag]);
    else
        n1 = t1.nstates;
        n2 = t2.nstates;
        fout = @(varargin)t1.fHandle(varargin{1:n1}).*t2.fHandle(varargin{n1+1:end});
        tout = term('fHandle',fout,'nstates',n1+n2);
    end
end