
function [W,H,runtime,niter] = nmf_wrap(v, model, verbose)
if nargin<3,
 verbose=false;
end
[W,H,elapsed,user,sys,niter] = brunet(v, size(model.W,2), verbose, model.W, model.H);
runtime.user = user;
runtime.sys = sys;
runtime.elapsed = elapsed;
end
