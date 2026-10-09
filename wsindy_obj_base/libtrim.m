syms q

pt = mean(cell2mat(arrayfun(@(j)arrayfun(@(i)mean(Uobj{j}.Uobs{i}),1:nstates),(1:ntraj)','uni',0)),1);
ord = 5;

args = num2cell(pt);
tays = zeros(nstates*ord,length(tags));
for j=1:length(tags)
    tt = [];
    for i=1:nstates
        arg_temp = args;
        arg_temp{i} = q;
        t1 = sym2poly(taylor(lib.terms{j}.fHandle(arg_temp{:}),q,0,'order',ord))';
        if length(t1)<ord
            t1 = [zeros(ord-length(t1),1);t1];
        end
        tt = [tt;flipud(t1)];
    end
    tays(:,j) = tt;
end

%%

sp = tays~=0;
Xsparse = unique(sp','rows')';
ii = cell(size(Xsparse,2),1);
As = {};
for j=1:size(Xsparse,2)
    ii = find(~any(sp-Xsparse(:,j),1));
    if length(ii)>1
        G = tays(:,ii);
        for l=1:nstates
            k = find(Xsparse((l-1)*ord+(1:ord),j),1);
            if ~isempty(k)
                G((l-1)*ord+(1:ord),:) = G((l-1)*ord+(1:ord),:)./G((l-1)*ord+k,:);
            end
        end
        A = real(sqrt((vecnorm(G)').^2 - 2*(G')*G + (vecnorm(G).^2)));
        A = A-triu(A);
        A = and(A>0,A<0.5);
        if any(A(:))
            As = [As,{{Xsparse(:,j),ii,G,A}}];
        end
    end
end

%%

rm_inds = [];
for i=1:length(As)
    while any(As{i}{end}(:))
        r={[],[]};
        [r{:}] = find(As{i}{end});
        r = [r{:}];
        r = unique(r(:));
        [~,a] = max(sum(abs(As{i}{3}([ord:ord:end],r))));
        rm = r(a);
        rm_inds = [rm_inds As{i}{2}(rm)];
        As{i}{2} = As{i}{2}([1:rm-1 rm+1:end]);
        As{i}{3} = As{i}{3}(:,[1:rm-1 rm+1:end]);
        As{i}{4} = As{i}{4}(:,[1:rm-1 rm+1:end]);
        As{i}{4} = As{i}{4}([1:rm-1 rm+1:end],:);
    end
end

% J = length(lib.terms);
% lib.terms = lib.terms(setdiff(1:J,rm_inds)); 
% lib.tags = lib.tags(setdiff(1:J,rm_inds));
% 










