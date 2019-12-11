function [bestPath,maxPathProb,T1,T2] = myViterbi(spikes,PI,A,B,dt)
if (size(A,1) ~= size(A,2))
    error('Transition Matrix is not square')
else
    nStates = size(A,1);
end
poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt);
nTimeSteps = size(spikes,2);

T1 = zeros(nStates,nTimeSteps);
T2 = zeros(nStates,nTimeSteps);
for i=1:nStates
    T1(i,1) = log(PI(i)) + log(prod(poiss(B(:,i),spikes(:,1))));
    T2(i,1) = 0;
end

for j=2:nTimeSteps
    for i=1:nStates
        vec1 = T1(:,j-1) + log(A(:,i)) + log(prod(poiss(B(:,i),spikes(:,j))));
        vec2 = T1(:,j-1) + log(A(:,i));
        T1(i,j) = max(vec1);
        [~,ind] = max(vec2);
        T2(i,j) = ind;
    end
end

[maxPathProb,bestPathEndState] = max(T1(:,end));
bestPath = zeros(1,nTimeSteps);
bestPath(end) = bestPathEndState;
for t=fliplr(1:(nTimeSteps-1))
    bestPath(t) = T2(bestPath(t+1),t+1);
end
end

