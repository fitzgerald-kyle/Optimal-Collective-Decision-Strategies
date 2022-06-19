function [X,over,dT] = diffusionTrial_oneAgent(p, agent) 
% p contains parameters. over is true if the agent decided on (+)).
    dt = p.dt; D = p.D;
    if agent == 1
        mu = p.mu1; H = p.H1; L = p.L1;  
    else
        mu = p.mu2; H = p.H2; L = p.L2;
    end
    
    X = zeros(2,1); % decision states

    i = 2;
    over = 0; under = 0;
    while ~over && ~under
        X(i) = X(i-1) + mu*dt + sqrt(2*D)*sqrt(dt)*randn;
        over = X(i) >= H; under = X(i) <= L;
        i = i+1;
    end
    
    dT = dt*(length(X)-1);
end