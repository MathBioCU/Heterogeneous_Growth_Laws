function out=arrayfunvec(arr_in,f,d,uni)
    if ~exist('uni','var')
        uni = 1;
    end
    dims=size(arr_in);
    dim = length(dims);
    e = ones(dim);
    e(1:dim+1:end)=dims;
    e=num2cell(e);
    inds=arrayfun(@(x)ones(e{x:dim:end}),1:dim,'uni',0);
    inds{d} = dims(d);
    arr_in_cell = mat2cell(arr_in,inds{:});
    out = cellfun(@(x)f(x),arr_in_cell,'uni',0);
    if and(all(cellfun(@(A)isequal(size(out{1}),size(A)),out(2:end))),uni)
        out = cell2mat(out);
    end
end