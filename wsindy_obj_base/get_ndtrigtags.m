% f0 = [];
% fmax = 5;
% nstates = 3;
% [tags,s] = get_ndtrigtags(nstates,fmax,f0);
% plot3(s(:,1),s(:,2),s(:,3),'o')

function [s,tags] = get_ndtrigtags(nstates,fmax,f0,nrm,getfuns)

    if isempty(f0)
        f0 = ones(1,nstates);
    elseif length(f0)<nstates
        f0 = ones(1,nstates)*f0(1);
    end
    
    x = repmat({-fmax:fmax},1,nstates-1);
    x = [x,{x{1}(1:fmax+1)}];
    s = cell(1,nstates);
    [s{:}] = ndgrid(x{:});
    s = cell2mat(cellfun(@(x)x(:),s,'uni',0));
    if isequal(nrm,inf)
        s = unique(s,'rows');
    elseif isequal(nrm,1)
        s = unique(s(sum(abs(s),2)<=fmax,:),'rows');
    end
    s = unique_to_sign_array(s);
    s = s.*f0;

    tags = {};
    if getfuns==1
        vars = arrayfun(@(i)['x',num2str(i)],1:nstates,'uni',0);
        expr = [];
        for i=1:nstates
            expr = [expr,vars{i}];
            if i < nstates
                expr = [expr,','];
            end
        end
        for j=1:size(s,1)
            str1 = ['@(',expr,')cos(',num2str(s(j,1)),'*',vars{1}];
            str2 = ['@(',expr,')sin(',num2str(s(j,1)),'*',vars{1}];
            for k=2:nstates
                if s(j,k)
                    if s(j,k)<0
                        if s(j,k) == -1
                            str1 = [str1,'-',vars{k}];
                            str2 = [str2,'-',vars{k}];
                        else
                            str1 = [str1,num2str(s(j,k)),'*',vars{k}];
                            str2 = [str2,num2str(s(j,k)),'*',vars{k}];
                        end
                    else
                        if s(j,k)==1
                            str1 = [str1,'+',vars{k}];
                            str2 = [str2,'+',vars{k}];
                        else
                            str1 = [str1,'+',num2str(s(j,k)),'*',vars{k}];
                            str2 = [str2,'+',num2str(s(j,k)),'*',vars{k}];
                        end
                    end
                end
            end
            str1 = [str1,')'];
            str2 = [str2,')'];
            tags = [tags;{eval(str1)}];
            if any(s(j,:))
                tags = [tags;{eval(str2)}];
            end
        end
    end

    s1 = s*(1+1i);
    ll = find(all(s==0,2));
    if isempty(ll)
        s2 = s*(1-1i);
    else
        s2 = s([1:ll-1 ll+1:end],:)*(1-1i);
    end
    s = [s1;s2];

end

function s=unique_to_sign_array(s)
    i=1;
    while i<size(s,1)
        j=i+1;
        check = 1;
        while and(check,j<=size(s,1))
            check = unique_to_sign(s(i,:),s(j,:));
            j=j+1;
        end
        if ~check
            s = s([1:j-2 j:end],:);
        end
        i=i+1;
    end
end

function bool=unique_to_sign(k,j)
    if or(all(k==j),all(k==-j))
        bool=false;
    else
        bool=true;
    end
end

