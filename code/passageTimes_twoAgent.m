function [PT1H, PT1L, PT2H, PT2L] = passageTimes_twoAgent(p, N, simple, varargin)
% Returns first passage times for N diffusion trials for the two-agent
% scenario.
%
% varargin contains the identity of the agent on which first passage is
% conditioned.

    PT1H = []; % upper boundary passage times for agent 1, etc
    PT1L = [];
    PT2H = [];
    PT2L = [];
    for i=1:N
        if simple
            [X1,X2,over1,over2,firstAgent] = diffusionTrial_twoAgent_simple(p);
        else
            [X1,X2,over1,over2,firstAgent] = diffusionTrial_twoAgent(p,1);
        end
        if nargin == 3 || (varargin{1} == 1 && firstAgent == 1)
            if over1
                PT1H(end+1) = p.dt*(length(X1)-1);
            else
                PT1L(end+1) = p.dt*(length(X1)-1);
            end
        end
        if nargin == 3 || (varargin{1} == 2 && firstAgent == 2)
            if over2
                PT2H(end+1) = p.dt*(length(X2)-1);
            else
                PT2L(end+1) = p.dt*(length(X2)-1);
            end
        end
    end
    
    PT1H = sort(PT1H); PT1L = sort(PT1L);
    PT2H = sort(PT2H); PT2L = sort(PT2L);
end