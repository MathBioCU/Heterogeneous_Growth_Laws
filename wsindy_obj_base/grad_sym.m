function grads = grad_sym(features,nstates)
    J = length(features);
    grads = cell(J,1);
    if J>0
        if ~exist('nstates','var')
            a = strfind(functions(features{1}).function,')');
            a = a(1);
            nstates = 1+length(strfind(functions(features{1}).function(1:a),','));
        end
        X = str2sym(strcat('x',num2str((1:nstates)')));
        Xc = sym2cell(X);
        for j=1:J
            f = features{j};
            g = gradient(f(Xc{:}),X);
            grads{j} = matlabFunction(g,'vars',Xc);
        end
    end
end