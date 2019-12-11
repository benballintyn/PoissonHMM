function [beta] = backward(spikes,nStates,dt,A,B,norms)
poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt);
nTimeSteps = size(spikes,2);
beta = zeros(nStates,nTimeSteps);
beta(:,end) = 1;
for t=fliplr(1:(nTimeSteps-1))
    for s=1:nStates
        beta(s,t) = sum(beta(:,t+1).*A(s,:)'.*prod(poiss(B,spikes(:,t+1)))');
    end
    beta(:,t) = beta(:,t)./norms(t);
end
end

