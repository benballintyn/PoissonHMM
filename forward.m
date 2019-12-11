function [alpha,norms] = forward(spikes,nStates,dt,PI,A,B)
poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt);
nTimeSteps = size(spikes,2);
for i=1:nStates
    alpha(i,1) = PI(i)*prod(poiss(B(:,i),spikes(:,1)));
end
norms(1) = sum(alpha(:,1));
alpha(:,1) = alpha(:,1)./norms(1);
for t=2:nTimeSteps
    for s=1:nStates
        alpha(s,t) = prod(poiss(B(:,s),spikes(:,t)))*sum(alpha(:,t-1).*A(:,s));
    end
    norms(t) = sum(alpha(:,t));
    alpha(:,t) = alpha(:,t)./norms(t);
end
end

