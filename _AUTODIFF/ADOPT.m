function [x,m,v,nIter] = ADOPT(grad,x,beta1,beta2,step_adam,epsilon,m,v,nIter)
warning('ADOPT implementation might be wrong please double check. did not use the clamp that did not work. needs a bigger epsilon to converge than adam but beta2 value is free')
%beta1 = 0.9;
%beta2 = 0.999;
%epsilon = eps(single(1));
if nIter==1
    v=grad.^2;
    m=0;
end
%nIter.^0.25
%m=beta1*m+(1-beta1).*clip(grad./max(sqrt(v),epsilon),-nIter.^0.25,-nIter.^0.25);
m=beta1*m+(1-beta1).*grad./max(sqrt(v),epsilon);
x = x + step_adam*m;
v = beta2.*v + (1 - beta2) .* (grad.^2);
nIter=nIter+1;

%{
m = beta1.*m + (1 - beta1) .* grad;
v = beta2.*v + (1 - beta2) .* (grad.^2);

mHat = m./(1 - beta1^nIter);
vHat = v./(1 - beta2^nIter);

x = x + step_adam*mHat./(sqrt(vHat) + epsilon);
nIter=nIter+1;
%}
end

