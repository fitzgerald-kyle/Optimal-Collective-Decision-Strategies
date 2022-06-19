function [X1,X2,over1,over2,instant] = diffusionTrial_twoAgent(p,firstStop) 
% firstStop is a boolean indicating whether to stop computing after the 
% first decision is made (to save computing time).

% over[1/2] is true if agent [1/2] decided at its upper threshold. instant 
% is a boolean indicating whether the second decision was made 
% instantaneously after the first. 

    dt = p.dt; mu1 = p.mu1; mu2 = p.mu2; D = p.D; qp = p.qp; qn = p.qn;
    H1 = p.H1; L1 = p.L1; H2 = p.H2; L2 = p.L2; 

    X1 = zeros(2,1); X2 = zeros(2,1); % decision states

    i = 2;
    over1 = 0; under1 = 0; over2 = 0; under2 = 0; instant = 0;
    while ~over1 && ~under1 && ~over2 && ~under2
        X1(i) = X1(i-1) + mu1*dt + sqrt(2*D)*sqrt(dt)*randn;
        X2(i) = X2(i-1) + mu2*dt + sqrt(2*D)*sqrt(dt)*randn;
        over1 = X1(i) >= H1; under1 = X1(i) <= L1;
        over2 = X2(i) >= H2; under2 = X2(i) <= L2;
        i = i+1;
    end
    
    if over1 || under1
        if firstStop
            return
        end
        X2(i-1) = X2(i-1) + qp*over1 - qn*under1;
        over2 = X2(i-1) >= H2; under2 = X2(i-1) <= L2;
        instant = over2 || under2;
        while ~over2 && ~under2
            X2(i) = X2(i-1) + mu2*dt + sqrt(2*D)*sqrt(dt)*randn;
            over2 = X2(i) >= H2; under2 = X2(i) <= L2;
            i = i+1;
        end
    else
        if firstStop
            return
        end
        X1(i-1) = X1(i-1) + qp*over2 - qn*under2;
        over1 = X1(i-1) >= H1; under1 = X1(i-1) <= L1;
        instant = over1 || under1;
        while ~over1 && ~under1
            X1(i) = X1(i-1) + mu1*p.dt + sqrt(2*D)*sqrt(dt)*randn;
            over1 = X1(i) >= H1; under1 = X1(i) <= L1;
            i = i+1;
        end
    end
end