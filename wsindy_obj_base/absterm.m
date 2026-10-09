classdef absterm < handle
    properties
        fHandle
        ftag
        gradterms
        gradon
        nstates
        linOp
    end
    
    methods
       evalterm(obj)
       evalgrads(obj)
       diffmat(obj)
    end

end