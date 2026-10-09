function tnew=composeterm(told,h,hp)
    tnew = term('fHandle',@(varargin) h(told.fHandle(varargin{:})));
    tnew.gradterms = cell(1,told.numargs);
    tnew.numargs = told.numargs;
    for j=1:t1.numargs
        tnew.gradterms{j} = term('fHandle',@(varargin) hp(told.fHandle(varargin{:})).*told.gradterms{j}.fHandle(varargin{:}));
    end
    tnew.tag = {told.tag,h}; %%% need better way to combine these
end