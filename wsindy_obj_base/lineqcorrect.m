function w = lineqcorrect(G,b,A,d)
    w = G \ b;
    if ~any([isempty(A) isempty(b)])
        % H = (G'*G) \ A';
        H = lsqminnorm(G'*G, A');
        e = A*H;
        if ~isequal(e,0)
            mu = e \ (A*w-d);
            w = w - H*mu;
        end
    end
end