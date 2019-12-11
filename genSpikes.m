function [spikes] = genSpikes(rates,stateSeq,dt)
spikes = zeros(size(rates,1),length(stateSeq)/dt);
for i=1:size(spikes,1)
    for j=1:length(stateSeq)
        curRate = rates(i,stateSeq(j));
        v = rand(1,(1/dt)); v = v < curRate*dt;
        spikes(i,((j-1)*(1/dt)+1):((j-1)*(1/dt) + 1/dt)) = v;
    end
end
end

