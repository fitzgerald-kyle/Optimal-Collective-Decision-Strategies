function [PTH, PTL] = passageTimes_oneAgent(p, N, agent)
% Returns first passage times for N diffusion trials for the one-agent
% scenario.

    PTH = [];
    PTL = [];
    for i=1:N
        [~,over,dT] = diffusionTrial_oneAgent(p, agent);
        
        if over
            PTH(end+1) = dT;
        else
            PTL(end+1) = dT;
        end
    end
    
    PTH = sort(PTH); PTL = sort(PTL);
end