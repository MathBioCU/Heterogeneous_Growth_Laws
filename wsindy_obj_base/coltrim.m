function inds = coltrim(A,lam,y,exinds)
    if ~exist('exinds','var')
        exinds = [];
    end
    if isequal(exinds,'all')
        exinds = 1:size(A,2);
    end
    B = abs(((A'*A-diag(diag(A'*A)))./vecnorm(A))./(vecnorm(A)'));
    C = B > lam;
    proj = abs(y'*A./vecnorm(A));
    inds = true(size(B,1),1);
   
    C(vecnorm(A)==0,:) = 0;
    C(:,vecnorm(A)==0) = 0;
    inds(vecnorm(A)==0) = false;

    meth = 2;
    k=1;
    while all([any(C(:)), k<size(B,1)*(size(B,1)-1), length(find(inds))>2])
    
        if meth == 0 % take out minimum projection term from set of all correlated terms
            inds_corr = find(sum(C,2)~=0);
            [~,a] = min(proj(inds_corr));
            a = inds_corr(a);
            
        elseif meth == 1 % take out term with maximum sum correlations
            [~,a] = max(sum(C,2));
    
        elseif meth == 2 % find maximally correlated pair, take out term with lower projection
            [r, ~] = find(B == max(B(C)));
            [~,a] = min(proj(r));
            if ismember(r(a),exinds)
                a = r(max(mod(a+1,2),1));
            else
                a = r(a);
            end            
        end
    
        if ~ismember(a,exinds)
            C(:,a) = 0;
            C(a,:) = 0;
            inds(a) = false;
        end

        % imagesc(C)
        % drawnow
        k=k+1;
    end
    
    inds = find(inds);

end