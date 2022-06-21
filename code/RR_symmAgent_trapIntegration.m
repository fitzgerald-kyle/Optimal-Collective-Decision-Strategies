function RR = RR_symmAgent_trapIntegration(p, X, TI, N, ints)
% (Uses the trapezoidal rule to compute reward rate.)

% Calculates expected reward rate for two agents in consecutive
% environments with random probability of positive or negative associated
% drift, where decisions are separated by time TI. Agents are "symmetric"
% in the sense that they set decision thresholds at the same (possibly
% asymmetric) values, they receive the same (possibly asymmetric) "kick" 
% from the other agent when it is the first to decide at a boundary, and
% they receive the same (possibly asymmetric) reward upon decision at the
% boundary that matches the drift of the current environment. For 
% simplicity, we assume that agents are in a mu=1 or mu=-1 environment each
% with probability 1/2. Due to agent symmetry, WLOG we perform all
% calculations for one agent, conditioning on the other agent deciding 
% first if applicable. N is the number of image terms to use in the
% probability densities that solve the Focker-Planck equation. ints is 
% the number of trapezoidal integration intervals.

% set parameters contained in 'X' explicitly for maximization routine
p = parameters('p', p, 'H1', X(1), 'L1', X(2), 'qp1', X(3), 'qn1', X(4)); 
r = exp_reward(p, N, ints); 
TL = exp_TL(p, N, ints);
RR = 2*r / (TL + TI); 
end


function r = exp_reward(p, N, ints)
% Expected reward for each agent across all environments, WLOG agent 1.

probAtUpper = 1/2 * ( P_crossBefore(p, p.H1, 1, N, ints) + ...
    P_instCrossAfter(p, p.H1, 1, N, ints) + ...
    P_diffCrossAfter(p, p.H1, 1, N, ints) );
probAtLower = 1/2 * ( P_crossBefore(p, p.L1, -1, N, ints) + ...
    P_instCrossAfter(p, p.L1, -1, N, ints) + ...
    P_diffCrossAfter(p, p.L1, -1, N, ints) );

r = p.R1p/2 * probAtUpper + p.R1n/2 * probAtLower;
end


function TL = exp_TL(p, N, ints)
% Expected time for the last agent to make a decision (WLOG for agent 2 to
% decide, conditioned on agent 1 deciding first).

ti = linspace(0,1,ints+1); hi = 1/ints;
denomPos1 = 0; denomNeg1 = 0; numPos1 = 0; numNeg1 = 0;
for i= 2:ints % don't go to ints+1 because last term DNE and is zero anyway
    tiloc = ti(i)/(1-ti(i));
    tj = linspace(0,tiloc,ints+1); hj = tiloc/ints;
    sums_pos1 = (-p.D*dcdx(p,p.H1,tiloc,1,N,'mu',1) + ...
        p.D*dcdx(p,p.L1,tiloc,1,N,'mu',1))/2;
    sums_neg1 = (-p.D*dcdx(p,p.H1,tiloc,1,N,'mu',-1) + ...
        p.D*dcdx(p,p.L1,tiloc,1,N,'mu',-1))/2;
    for j= 2:ints
        sums_pos1 = sums_pos1 + ...
            (-p.D*dcdx(p,p.H1,tj(j),1,N,'mu',1) + p.D*dcdx(p,p.L1,tj(j),1,N,'mu',1));
        sums_neg1 = sums_neg1 + ...
            (-p.D*dcdx(p,p.H1,tj(j),1,N,'mu',-1) + p.D*dcdx(p,p.L1,tj(j),1,N,'mu',-1));
    end
    sums_pos1 = sums_pos1*hj;
    sums_neg1 = sums_neg1*hj;
    numPos1 = numPos1 - tiloc*(p.D*dcdx(p,p.H1,tiloc,1,N,'mu',1) ...
        + p.D*dcdx(p,p.L1,tiloc,1,N,'mu',1)) * sums_pos1 / (1-ti(i))^2;
    numNeg1 = numNeg1 - tiloc*(p.D*dcdx(p,p.H1,tiloc,1,N,'mu',1) ...
        + p.D*dcdx(p,p.L1,tiloc,1,N,'mu',1)) * sums_neg1 / (1-ti(i))^2;
    denomPos1 = denomPos1 - (p.D*dcdx(p,p.H1,tiloc,1,N,'mu',1) ...
        + p.D*dcdx(p,p.L1,tiloc,1,N,'mu',1)) * sums_pos1 / (1-ti(i))^2;
    denomNeg1 = denomNeg1 - (p.D*dcdx(p,p.H1,tiloc,1,N,'mu',1) ...
        + p.D*dcdx(p,p.L1,tiloc,1,N,'mu',1)) * sums_neg1 / (1-ti(i))^2;
end
numPos1 = numPos1*hi; numNeg1 = numNeg1*hi;
denomPos1 = denomPos1*hi; denomNeg1 = denomNeg1*hi;

TL = 1/2*(numPos1/denomPos1 + numNeg1/denomNeg1);
end

function P = P_crossBefore(p, thresh, mu, N, ints)
% Probability that one agent crosses the decision threshold at "thresh",
% given that it decides before the other agent.

ti = linspace(0,1,ints+1); hti = 1/ints;
tj = linspace(0,1,ints+1); htj = 1/ints;
P = 0;
for i= 2:ints % don't go to ints+1 because last term DNE and is zero anyway
    tiloc = ti(i)/(1-ti(i));
    sumstj = (-p.D*dcdx(p,p.H1,tiloc,1,N,'mu',mu) + ...
            p.D*dcdx(p,p.L1,tiloc,1,N,'mu',mu))/2;
    for j= 2:ints
        tjloc = tj(j)/(1-tj(j));
        sumstj = sumstj + (-p.D*dcdx(p,p.H1,tiloc+tjloc,1,N,'mu',mu) + ...
            p.D*dcdx(p,p.L1,tiloc+tjloc,1,N,'mu',mu)) / (1-tj(j))^2;
    end
    sumstj = sumstj*htj;
    P = P - sign(thresh)*p.D*dcdx(p,thresh,tiloc,1,N,'mu',mu) * ...
        sumstj / (1-ti(i))^2;
end
P = P*hti;
end

function P = P_instCrossAfter(p, thresh, mu, N, ints)
% Probability that one agent is instantaneously "kicked" across the 
% decision threshold at "thresh", given the other agent's crossing of 
% "thresh" before.

ti = linspace(0,1,ints+1); hti = 1/ints;
tj = linspace(0,1,ints+1); htj = 1/ints;
if sign(thresh) == 1
    xjnum = linspace(p.H1-p.qp1,p.H1,ints+1); hxjnum = p.qp1/ints;
else
    xjnum = linspace(p.L1,p.L1+p.qn1,ints+1); hxjnum = p.qn1/ints;
end
xjLtoH = linspace(p.L1,p.H1,ints+1); hxjLtoH = (p.H1-p.L1)/ints;
P = 0;
for i= 2:ints % don't go to ints+1 because last term DNE and is zero anyway
    tiloc = ti(i)/(1-ti(i));
    sumstj = (-p.D*dcdx(p,p.H1,tiloc,1,N,'mu',mu) + ...
        p.D*dcdx(p,p.L1,tiloc,1,N,'mu',mu))/2; % first term
    sumsxjnum = (c(p,xjnum(1),tiloc,1,N,'mu',mu) + ...
        c(p,xjnum(end),tiloc,1,N,'mu',mu))/2; % first and last terms
    sumsxjden = (c(p,xjLtoH(1),tiloc,1,N,'mu',mu) + ...
        c(p,xjLtoH(end),tiloc,1,N,'mu',mu))/2; % first and last terms
    for j= 2:ints
        tjloc = tj(j)/(1-tj(j));
        sumstj = sumstj + (-p.D*dcdx(p,p.H1,tiloc+tjloc,1,N,'mu',mu) + ...
            p.D*dcdx(p,p.L1,tiloc+tjloc,1,N,'mu',mu)) / (1-tj(j))^2;
        sumsxjnum = sumsxjnum + c(p,xjnum(j),tiloc,1,N,'mu',mu);
        sumsxjden = sumsxjden + c(p,xjLtoH(j),tiloc,1,N,'mu',mu);
    end
    sumstj = sumstj*htj;
    sumsxjnum = sumsxjnum*hxjnum; sumsxjden = sumsxjden*hxjLtoH;
    P = P - sign(thresh)*p.D*dcdx(p,thresh,tiloc,1,N,'mu',mu) * ...
        sumstj * sumsxjnum / sumsxjden / (1-ti(i))^2;
end
P = P*hti;
end

function P = P_diffCrossAfter(p, thresh, mu, N, ints)
% Probability that one agent diffuses across the decision threshold at 
% "thresh", given the other agent's crossing of one of the thresholds before.

if sign(thresh) == 1
    epsFunc = @(x) (exp(-mu*x/p.D)-exp(-mu*p.L1/p.D)) / ...
        (exp(-mu*p.H1/p.D)-exp(-mu*p.L1/p.D));
else 
    epsFunc = @(x) (exp(-mu*p.H1/p.D)-exp(-mu*x/p.D)) / ...
        (exp(-mu*p.H1/p.D)-exp(-mu*p.L1/p.D));
end

ti = linspace(0,1,ints+1); hti = 1/ints;
tj = linspace(0,1,ints+1); htj = 1/ints;
xjnum1 = linspace(p.L1,p.H1-p.qp1,ints+1); hxjnum1 = (p.H1-p.qp1-p.L1)/ints;
xjnum2 = linspace(p.L1+p.qn1,p.H1,ints+1); hxjnum2 = (p.H1-p.qn1-p.L1)/ints;
xjLtoH = linspace(p.L1,p.H1,ints+1); hxjLtoH = (p.H1-p.L1)/ints;
P = 0;
for i= 2:ints % don't go to ints+1 because last term DNE and is zero anyway
    tiloc = ti(i)/(1-ti(i));
    sumstj = (-p.D*dcdx(p,p.H1,tiloc,1,N,'mu',mu) + ...
        p.D*dcdx(p,p.L1,tiloc,1,N,'mu',mu))/2; % first term
    sumsxjnum1 = (c(p,xjnum1(1),tiloc,1,N,'mu',mu) + ...
        c(p,xjnum1(end),tiloc,1,N,'mu',mu))/2; % first and last terms
    sumsxjnum2 = (c(p,xjnum2(1),tiloc,1,N,'mu',mu) + ...
        c(p,xjnum2(end),tiloc,1,N,'mu',mu))/2; % first and last terms
    sumsxjnumEscape = (c(p,xjLtoH(1),tiloc,1,N,'mu',mu)*epsFunc(xjLtoH(1)) + ...
        c(p,xjLtoH(end),tiloc,1,N,'mu',mu)*epsFunc(xjLtoH(end)))/2; % first and last terms
    sumsxjden = (c(p,xjLtoH(1),tiloc,1,N,'mu',mu) + ...
        c(p,xjLtoH(end),tiloc,1,N,'mu',mu))/2; % first and last terms
    for j= 2:ints
        tjloc = tj(j)/(1-tj(j));
        sumstj = sumstj + (-p.D*dcdx(p,p.H1,tiloc+tjloc,1,N,'mu',mu) + ...
            p.D*dcdx(p,p.L1,tiloc+tjloc,1,N,'mu',mu)) / (1-tj(j))^2;
        sumsxjnum1 = sumsxjnum1 + c(p,xjnum1(j),tiloc,1,N,'mu',mu);
        sumsxjnum2 = sumsxjnum2 + c(p,xjnum2(j),tiloc,1,N,'mu',mu);
        sumsxjnumEscape = sumsxjnumEscape + ...
            c(p,xjLtoH(j),tiloc,1,N,'mu',mu)*epsFunc(xjLtoH(j));
        sumsxjden = sumsxjden + c(p,xjLtoH(j),tiloc,1,N,'mu',mu);
    end
    sumstj = sumstj*htj;
    sumsxjnum1 = sumsxjnum1*hxjnum1; sumsxjnum2 = sumsxjnum2*hxjnum2;
    sumsxjnumEscape = sumsxjnumEscape*hxjLtoH;
    sumsxjden = sumsxjden*hxjLtoH;
    P = P + (-p.D*dcdx(p,p.H1,tiloc,1,N,'mu',mu)*sumsxjnum1 + ...
        p.D*dcdx(p,p.L1,tiloc,1,N,'mu',mu)*sumsxjnum2) ...
        * sumstj * sumsxjnumEscape / sumsxjden^2  / (1-ti(i))^2;
end
P = P*hti;
end