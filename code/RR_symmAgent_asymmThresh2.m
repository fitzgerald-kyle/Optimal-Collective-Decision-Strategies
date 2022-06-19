function RR = RR_symmAgent_asymmThresh2(p, X, TI, N)
% Calculates expected reward rate for two agents in consecutive
% environments with random probability of positive or negative associated
% drift, where decisions are separated by time TI. Agents are "symmetric"
% in the sense that they set decision thresholds at the same (possibly
% asymmetric) values, they receive the same "kick" from the other 
% agent when it is the first to decide at a given boundary, and
% they receive the same (possibly asymmetric) reward upon decision at the
% boundary that matches the drift of the current environment. For 
% simplicity, we assume that agents are in a mu=1 or mu=-1 environment each
% with probability 1/2. Due to agent symmetry, WLOG we perform all
% calculations for one agent, conditioning on the other agent deciding 
% first if applicable. N is the number of image terms to use in the
% probability densities that solve the Focker-Planck equation.

% set parameters contained in 'X' explicitly for maximization routine
p = parameters('p', p, 'H1', X(1), 'L1', X(2), 'qp1', X(3), 'qn1', X(4));

%tic;
r = exp_reward(p, N);
%T=toc;
TL = exp_TL(p, N);
RR = 2*r / (TL + TI); 
end

function r = exp_reward(p, N)
% Expected reward for each agent across all environments, WLOG agent 1.
fH_cond_inf_P = @(t1,t2) fH(t1,1,p,N).*(fH(t2,1,p,N)+fL(t2,1,p,N));
fL_cond_inf_P = @(t1,t2) fL(t1,1,p,N).*(fH(t2,1,p,N)+fL(t2,1,p,N));
fH_cond_inf_N = @(t1,t2) fH(t1,-1,p,N).*(fH(t2,-1,p,N)+fL(t2,-1,p,N));
fL_cond_inf_N = @(t1,t2) fL(t1,-1,p,N).*(fH(t2,-1,p,N)+fL(t2,-1,p,N));
epsH_P = @(x) (exp(-x.*p.D)-exp(-p.L1.*p.D)) ./ (exp(-p.H1.*p.D)-exp(-p.L1.*p.D));
epsL_N = @(x) (exp(p.H1.*p.D)-exp(x.*p.D)) ./ (exp(p.H1.*p.D)-exp(p.L1.*p.D));

probAtUpper = 1/2 * ( P_crossBefore(fH_cond_inf_P) + ...
    P_instCrossAfter(p, fH_cond_inf_P, p.H1, 1, N) + ...
    P_diffCrossAfter(p, fH_cond_inf_P, fL_cond_inf_P, epsH_P, 1, N) );
probAtLower = 1/2 * ( P_crossBefore(fL_cond_inf_N) + ...
    P_instCrossAfter(p, fL_cond_inf_N, p.L1, -1, N) + ...
    P_diffCrossAfter(p, fH_cond_inf_N, fL_cond_inf_N, epsL_N, -1, N) );

r = p.R1p/2 * probAtUpper + p.R1n/2 * probAtLower;
end

function TL = exp_TL(p, N)
% Expected time for the last agent to make a decision (WLOG for agent 2 to
% decide, conditioned on agent 1 deciding first).
%tic;
tf = 10;
TL = integral2(@(t1,t2) t1.*fH(t1,1,p,N).*fH(t2,1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
    integral2(@(t1,t2) t1.*fH(t1,1,p,N).*fL(t2,1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
    integral2(@(t1,t2) t1.*fL(t1,1,p,N).*fH(t2,1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
    integral2(@(t1,t2) t1.*fL(t1,1,p,N).*fL(t2,1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
    integral2(@(t1,t2) t1.*fH(t1,-1,p,N).*fH(t2,-1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
    integral2(@(t1,t2) t1.*fH(t1,-1,p,N).*fL(t2,-1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
    integral2(@(t1,t2) t1.*fL(t1,-1,p,N).*fH(t2,-1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
    integral2(@(t1,t2) t1.*fL(t1,-1,p,N).*fL(t2,-1,p,N),0,tf,0,@(t1)t1,'AbsTol',1e-6);
%TL = integral2(@(t1,t2) t1.*((fH(t1,1,p,N)+fL(t1,1,p,N)).*(fH(t2,1,p,N)+fL(t2,1,p,N)) + ...
%    (fH(t1,-1,p,N)+fL(t1,-1,p,N)).*(fH(t2,-1,p,N)+fL(t2,-1,p,N))), ...
%    0,tf,0,@(t1)t1,'AbsTol',1e-6);
%TL = integral2(@(t1,t2) t1.*(fH(t1,1,p,N)+fL(t1,1,p,N)).*(fH(t2,1,p,N)+fL(t2,1,p,N)), ...
%    0,tf,0,@(t1)t1,'AbsTol',1e-6) + ...
%    integral2(@(t1,t2) t1.*(fH(t1,-1,p,N)+fL(t1,-1,p,N)).*(fH(t2,-1,p,N)+fL(t2,-1,p,N)), ...
%    0,tf,0,@(t1)t1,'AbsTol',1e-6);
%T4 = toc
end

function P = P_crossBefore(fthresh_cond_inf)
% Probability that one agent crosses the decision threshold at "thresh",
% given that it decides before the other agent.
%tic;
tf = 10;%inf;
P = integral2(@(t1,t2) fthresh_cond_inf(t1,t2),0,tf,@(t1)t1,tf,'AbsTol',1e-6);
%T1 = toc
end

function P = P_instCrossAfter(p, fthresh_cond_inf, thresh, mu, N)
% Probability that one agent is instantaneously "kicked" across the 
% decision threshold at "thresh", given the other agent's crossing of 
% "thresh" before.
%tic;
tf = 10;%inf;
if sign(thresh) == 1
    P = integral2(@(t1,t2) fthresh_cond_inf(t1,t2).*(intc_x(p,thresh,t1,1,N,mu) ...
        - intc_x(p,thresh-p.qp1,t1,1,N,mu)) ./ rhoDenom(t1,mu,p,N), ...
        0,tf,@(t1)t1,tf,'AbsTol',1e-6);
elseif sign(thresh) == -1
    P = integral2(@(t1,t2) fthresh_cond_inf(t1,t2).*(intc_x(p,thresh+p.qn1,t1,1,N,mu) ...
        - intc_x(p,thresh,t1,1,N,mu)) ./ rhoDenom(t1,mu,p,N), ...
        0,tf,@(t1)t1,tf,'AbsTol',1e-6);
end
%T2 = toc
end

function P = P_diffCrossAfter(p, fH_cond_inf, fL_cond_inf, epsFun, mu, N)
% Probability that one agent diffuses across the decision threshold at 
% "thresh", given the other agent's crossing of one of the thresholds before.
%tic;
tf = 10;
P = integral2(@(t1,t2) fH_cond_inf(t1,t2).*(intc_x(p,p.H1-p.qp1,t1,1,N,mu) - ...
        intc_x(p,p.L1,t1,1,N,mu)) ./ rhoDenom(t1,mu,p,N).^2 .* ...
        integral(@(x) c(p,x,t1,1,N,mu).*epsFun(x),p.L1,p.H1,'ArrayValued',1), ...
        0,tf,@(t1)t1,tf,'AbsTol',1e-6) + ...
        integral2(@(t1,t2) fL_cond_inf(t1,t2).*(intc_x(p,p.H1,t1,1,N,mu) - ...
        intc_x(p,p.L1+p.qn1,t1,1,N,mu)) ./ rhoDenom(t1,mu,p,N).^2 .* ...
        integral(@(x) c(p,x,t1,1,N,mu).*epsFun(x),p.L1,p.H1,'ArrayValued',1), ...
        0,tf,@(t1)t1,tf,'AbsTol',1e-6);
%T3 = toc
end

function val = fH(t,mu,p,N)
    val = -p.D * dcdx(p,p.H1,t,1,N,mu);
end

function val = fL(t,mu,p,N)
    val = p.D * dcdx(p,p.L1,t,1,N,mu);
end

function val = rhoDenom(t,mu,p,N)
    val = intc_x(p,p.H1,t,1,N,mu) - intc_x(p,p.L1,t,1,N,mu);
    if abs(val)==0, val = eps; end
end