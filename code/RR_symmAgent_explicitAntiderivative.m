function RR = RR_symmAgent_explicitAntiderivative(p, X, TI, N)
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
% probability densities that solve the Focker-Planck equation.

% set parameters contained in 'X' explicitly for maximization routine
if length(X)==4
    p = parameters('p', p, 'H1', X(1), 'L1', X(2), 'qp1', X(3), 'qn1', X(4));
elseif length(X)==2
    p = parameters('p', p, 'H1', X(1), 'L1', X(2));
end

r = exp_reward(p, N);
TL = exp_TL(p, N);
RR = 2*r / (TL + TI); 
end

function r = exp_reward(p, N)
% Expected reward for each agent across all environments, WLOG agent 1.
fH_cond_inf_P = @(t) -fH(t,1,p,N).*(intf_t(p,p.H1,t,N,1)+intf_t(p,p.L1,t,N,1));
fL_cond_inf_P = @(t) -fL(t,1,p,N).*(intf_t(p,p.H1,t,N,1)+intf_t(p,p.L1,t,N,1));
fH_cond_inf_N = @(t) -fH(t,-1,p,N).*(intf_t(p,p.H1,t,N,-1)+intf_t(p,p.L1,t,N,-1));
fL_cond_inf_N = @(t) -fL(t,-1,p,N).*(intf_t(p,p.H1,t,N,-1)+intf_t(p,p.L1,t,N,-1));

probAtUpper = 1/2 * ( P_crossBefore(fH_cond_inf_P) + ...
    P_instCrossAfter(p, fH_cond_inf_P, p.H1, 1, N) + ...
    P_diffCrossAfter(p, fH_cond_inf_P, fL_cond_inf_P, p.H1, 1, N) );
probAtLower = 1/2 * ( P_crossBefore(fL_cond_inf_N) + ...
    P_instCrossAfter(p, fL_cond_inf_N, p.L1, -1, N) + ...
    P_diffCrossAfter(p, fH_cond_inf_N, fL_cond_inf_N, p.L1, -1, N) );

r = p.R1p/2 * probAtUpper + p.R1n/2 * probAtLower;
end

function TL = exp_TL(p, N)
% Expected time for the last agent to make a decision (WLOG for agent 2 to
% decide, conditioned on agent 1 deciding first).
ti = eps;
tf = inf;
TL = integral(@(t) t.*(fH(t,1,p,N)+fL(t,1,p,N)) .* (intf_t(p,p.H1,t,N,1) ...
    +intf_t(p,p.L1,t,N,1)-intf_t(p,p.H1,ti,N,1)-intf_t(p,p.L1,ti,N,1)),ti,tf,'AbsTol',1e-10) ...
     + integral(@(t) t.*(fH(t,-1,p,N)+fL(t,-1,p,N)) .* (intf_t(p,p.H1,t,N,-1) ...
    +intf_t(p,p.L1,t,N,-1)-intf_t(p,p.H1,ti,N,-1)-intf_t(p,p.L1,ti,N,-1)),ti,tf,'AbsTol',1e-10);
end

function P = P_crossBefore(fthresh_cond_inf)
% Probability that one agent crosses the decision threshold at "thresh",
% given that it decides before the other agent.
ti = eps;
tf = 100;
P = integral(@(t) fthresh_cond_inf(t),ti,tf,'AbsTol',1e-10);
end

function P = P_instCrossAfter(p, fthresh_cond_inf, thresh, mu, N)
% Probability that one agent is instantaneously "kicked" across the 
% decision threshold at "thresh", given the other agent's crossing of 
% "thresh" before.
ti = eps;
tf = inf;
if sign(thresh) == 1
    P = integral(@(t) fthresh_cond_inf(t).*(intc_x(p,thresh,t,N,mu) ...
        - intc_x(p,thresh-p.qp1,t,N,mu)) ./ rhoDenom(t,mu,p,N), ...
        ti,tf,'AbsTol',1e-10);
elseif sign(thresh) == -1
    P = integral(@(t) fthresh_cond_inf(t).*(intc_x(p,thresh+p.qn1,t,N,mu) ...
        - intc_x(p,thresh,t,N,mu)) ./ rhoDenom(t,mu,p,N), ...
        ti,tf,'AbsTol',1e-10);
end
end

function P = P_diffCrossAfter(p, fH_cond_inf, fL_cond_inf, thresh, mu, N)
% Probability that one agent diffuses across the decision threshold at 
% "thresh", given the other agent's crossing of one of the thresholds before.
ti = eps;
tf = 100;
if sign(thresh) > 0
    epsInt = @(t) intcEpsH_x(p,p.H1,t,N,mu)-intcEpsH_x(p,p.L1,t,N,mu);
elseif sign(thresh) < 0
    epsInt = @(t) intc_x(p,p.H1,t,N,mu) - intc_x(p,p.L1,t,N,mu) - ...
        intcEpsH_x(p,p.H1,t,N,mu) + intcEpsH_x(p,p.L1,t,N,mu);
end
P = integral(@(t) fH_cond_inf(t).*(intc_x(p,p.H1-p.qp1,t,N,mu) - ...
        intc_x(p,p.L1,t,N,mu)) ./ rhoDenom(t,mu,p,N).^2 .* epsInt(t), ...
        ti,tf,'AbsTol',1e-10) + ...
        integral(@(t) fL_cond_inf(t).*(intc_x(p,p.H1,t,N,mu) - ...
        intc_x(p,p.L1+p.qn1,t,N,mu)) ./ rhoDenom(t,mu,p,N).^2 .* epsInt(t), ...
        ti,tf,'AbsTol',1e-10);
end

function val = fH(t,mu,p,N)
    val = -p.D * dcdx(p,p.H1,t,N,mu);
end

function val = fL(t,mu,p,N)
    val = p.D * dcdx(p,p.L1,t,N,mu);
end

function val = rhoDenom(t,mu,p,N)
    val = intc_x(p,p.H1,t,N,mu) - intc_x(p,p.L1,t,N,mu);
    val(abs(val)<=sqrt(realmin)) = sqrt(realmin);
end

function cfunc = c(p,x,t,numTerms,mu)
    cfunc = 0;
    for n=1:numTerms
        L = abs(p.L1);
        len = L+p.H1;
        w = n*pi./len;
        k = mu^2/(4*p.D) + p.D*w.^2;
        A = 2*sin(n*pi*L./len) ./ len;
        cfunc = cfunc + A .* ...
            exp(-k.*t).*sin(w.*(x+L));
    end
    cfunc = cfunc .* exp(mu*x/(2*p.D));
end

function dcdxfunc = dcdx(p,x,t,numTerms,mu)
    dcdxfunc = 0;
    for n=1:numTerms
        L = abs(p.L1);
        len = L+p.H1;
        w = n*pi./len;
        k = mu^2/(4*p.D) + p.D*w.^2;
        A = 2*sin(n*pi*L./len) ./ len;
        dcdxfunc = dcdxfunc + A.*exp(-k.*t) .* ...
            ( mu*sin(w.*(x+L))/(2*p.D) + w.*cos(w.*(x+L)) );
    end
    dcdxfunc = dcdxfunc .* exp(mu*x/(2*p.D));
end

function intcfunc = intc_x(p,x,t,numTerms,mu)
    intcfunc = 0;
    for n=1:numTerms
        L = abs(p.L1);
        len = L+p.H1;
        w = n*pi./len;
        k = mu^2/(4*p.D) + p.D*w.^2;
        A = 2*sin(n*pi*L./len) ./ len;
        intcfunc = intcfunc + A.*exp(-k.*t) .* ...
            ( mu*sin(w.*(x+L)) - 2*p.D*w.*cos(w.*(x+L)) ) ./ ...
            (4*p.D^2*w.^2 + mu.^2);
    end
    intcfunc = intcfunc .* exp(mu*x/(2*p.D))*2*p.D;
end

function intf_func = intf_t(p,x,t,numTerms,mu)
    intf_func = 0;
    for n=1:numTerms
        L = abs(p.L1);
        len = L+p.H1;
        w = n*pi./len;
        k = mu^2/(4*p.D) + p.D*w.^2;
        A = 2*sin(n*pi*L./len) ./ len;
        intf_func = intf_func - A./k.*exp(-k.*t) .* ...
            ( mu*sin(w.*(x+L))./(2*p.D) + w.*cos(w.*(x+L)) );
    end
    intf_func = -intf_func .* exp(mu*x/(2*p.D)).*sign(x)*p.D;
end

function intcEpsHfunc = intcEpsH_x(p,x,t,numTerms,mu)
    intcEpsHfunc = 0;
    for n=1:numTerms
        L = abs(p.L1);
        len = L+p.H1;
        w = n*pi./len;
        k = mu^2/(4*p.D) + p.D*w.^2;
        A = 2*sin(n*pi*L./len) ./ len;
        intcEpsHfunc = intcEpsHfunc + A.*exp(-k.*t) .* ...
            ( 2*p.D*mu*(1+exp(mu*(x+L)/p.D)).*sin(w.*(x+L)) + ...
            4*p.D^2*w.*(1-exp(mu*(x+L)/p.D)).*cos(w.*(x+L)) ) ./ ...
            (4*p.D^2*w.^2 + mu.^2);
    end
    intcEpsHfunc = intcEpsHfunc .* exp(-mu*x/(2*p.D)) ...
        ./ ( exp(mu*L/p.D)-exp(-mu*p.H1/p.D) );
end