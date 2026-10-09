function I=gettrueinds(tags,true_tags)
    if isequal(class(tags),'double')
        tags = mat2cell(tags,ones(size(tags,1),1),size(tags,2));
    end
    if isequal(class(true_tags),'double')
        true_tags = mat2cell(true_tags,ones(size(true_tags,1),1),size(true_tags,2));
    end
    I = zeros(length(true_tags),1);
    for j=1:length(true_tags)
        ttemp = true_tags{j};
        if all(sum(ttemp~=0,1)<=1)
            ttemp = sum(ttemp,1);
        end
        for i=1:length(tags)
            tagtemp = tags{i};
            if all(sum(tagtemp~=0,1)<=1)
                tagtemp = sum(tagtemp,1);
            end
            if isequal(ttemp,tagtemp)
                I(j) = i;
            elseif all(size(ttemp)==size(tagtemp))
                try
                    if all(ismember(ttemp,tagtemp,'rows'))
                        I(j) = i;
                    end
                catch
                end
            elseif size(ttemp,1)==1
                [a,b] = ismember(tagtemp,ttemp,'rows');
                if and(any(a),all(tagtemp(~a,:)==0))
                    I(j) = i;
                end
            elseif size(tagtemp,1)==1
                [a,b] = ismember(ttemp,tagtemp,'rows');
                if and(any(a),all(ttemp(~a,:)==0))
                    I(j) = i;
                end
            else
                [a,b] = ismember(ttemp,tagtemp,'rows');
                if and(all(a),all(tagtemp(setdiff(1:size(tagtemp,1),b),:)==0))
                    I(j) = i;
                end
            end
        end
    end
end