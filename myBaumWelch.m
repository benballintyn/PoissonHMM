function [PI,A,B,alpha,beta,gamma,epsilons] = myBaumWelch(spikes,nStates,dt,maxIter)
nNeurons = size(spikes,1);
nTimeSteps = size(spikes,2);
poiss = @(lambda,n) ((lambda*dt).^n./factorial(n)).*exp(-lambda*dt);
minFR = (1/(nTimeSteps*dt));
% Initialize parameter estimates
PI = ones(nStates,1)*(1/nStates); % Initial state distribution
A = zeros(nStates,nStates); % Transition Matrix (columns = post, rows = pre)
for i=1:nStates
    for j=1:nStates
        if (i == j)
            A(i,j) = .99;
        else
            A(i,j) = .01/(nStates-1);
        end
    end
end

B = rand(nNeurons,nStates); % "Emission" (rate) matrix

notConverged = 1;
iterNum = 0;
while (notConverged && (iterNum < maxIter))
    % run forward algorithm
    [alpha,norms] = forward(spikes,nStates,dt,PI,A,B);
    % run backward algorithm
    [beta] = backward(spikes,nStates,dt,A,B,norms);
    % compute temporary variables
    gamma = zeros(nStates,nTimeSteps);
    epsilons = zeros(nStates,nStates,nTimeSteps-1);
    for t=1:nTimeSteps
        if (t < nTimeSteps)
            gamma(:,t) = (alpha(:,t).*beta(:,t))./sum(alpha(:,t).*beta(:,t));
            epsilonNumerator = zeros(nStates,nStates);
            for si=1:nStates
                for sj=1:nStates
                    epsilonNumerator(si,sj) = alpha(si,t)*A(si,sj)*beta(sj,t)*prod(poiss(B(:,sj),spikes(:,t+1)));
                end
            end
            epsilons(:,:,t) = epsilonNumerator./sum(epsilonNumerator(:));
        end
    end
    
    % stor old paramters for convergence check
    oldPI = PI;
    oldA = A;
    oldB = B;
    % update parameters
    PI = gamma(:,1);
    Anumer = sum(epsilons,3);
    Adenom = sum(gamma,2);
    A = Anumer./Adenom;
    A = A./sum(A,2);
    B = ((spikes*gamma')./sum(gamma,2)')/dt;
    B = max(minFR,B);
    % check for convergence
    notConverged = isConverged(oldPI,oldA,oldB,PI,A,B);
    iterNum=iterNum+1;
    disp(['done with iteration #' num2str(iterNum)])
end
end

function notConverged = isConverged(oldPI,oldA,oldB,PI,A,B)
dPI = sqrt(sum((oldPI - PI).^2));
dA = sqrt(sum((oldA(:) - A(:)).^2));
dB = sqrt(sum((oldB(:) - B(:)).^2));
disp(['dPI = ' num2str(dPI) ' dA = ' num2str(dA) ' dB = ' num2str(dB)])
thresh = 1e-4;
if (dPI < thresh && dA < thresh && dB < thresh)
    notConverged = 0;
else
    notConverged = 1;
end
end

