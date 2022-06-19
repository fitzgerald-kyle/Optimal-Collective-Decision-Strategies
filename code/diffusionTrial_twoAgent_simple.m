function [X1,X2,over1,over2,firstAgent] = diffusionTrial_twoAgent_simple(p) 
% over[1/2] is true if agent [1/2] decided at its upper threshold. 
% firstAgent is the identity of the first agent to make a decision.

    dt = p.dt; D = p.D; 
    mu1 = p.mu1; qp1 = p.qp1; qn1 = p.qn1; H1 = p.H1; L1 = p.L1;
    mu2 = p.mu2; qp2 = p.qp2; qn2 = p.qn2; H2 = p.H2; L2 = p.L2; 

    X1 = zeros(2,1); X2 = zeros(2,1); % decision states

    i = 2;
    over1 = 0; under1 = 0; over2 = 0; under2 = 0;
    while ~over1 && ~under1 && ~over2 && ~under2
        X1(i) = X1(i-1) + mu1*dt + sqrt(2*D)*sqrt(dt)*randn;
        X2(i) = X2(i-1) + mu2*dt + sqrt(2*D)*sqrt(dt)*randn;
        over1 = X1(i) >= H1; under1 = X1(i) <= L1;
        over2 = X2(i) >= H2; under2 = X2(i) <= L2;
        i = i+1;
    end
    
    
    if over1 || under1
        firstAgent = 1;
        X2(i-1) = X2(i-1) + qp2*over1 - qn2*under1;
        over2 = X2(i-1) >= (H2+L2)/2;
    else
        firstAgent = 2;
        X1(i-1) = X1(i-1) + qp1*over2 - qn1*under2;
        over1 = X1(i-1) >= (H1+L1)/2;
    end
end