function RR = RR_symmAgent_CrankNicolson_simple(p, X, TI)
% (Uses the Crank-Nicolson method to solve the Focker-Planck equation and
% determine an agent's expected reward rate.)

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
    p = parameters('p', p, 'H1', X(1), 'L1', X(2), 'qp1', X(3), 'qn1', X(4), ...
        'H2', X(1), 'L2', X(2), 'qp2', X(3), 'qn2', X(4));
elseif length(X)==2
    p = parameters('p', p, 'H1', X(1), 'L1', X(2));
end

frac = p.H1/(p.H1-p.L1);
numInt = 200;
rem = mod(numInt*frac, 1);
for i=201:300
    if rem==0
        break
    end
    if mod(i*frac, 1) < rem
        rem = mod(i*frac, 1);
        numInt = i;
    end
end

% Set grid spacing hx and ht, ensuring that the relationship between the
% two adheres to the CFL condition
hx = (p.H1-p.L1)/numInt;
min_char_t = 3*min([p.H1, abs(p.L1)]); 
tmax = 3*max([p.H1, abs(p.L1)]);
ht = min([min_char_t/1000, hx^2]);

x = (p.L1:hx:p.H1)';
t = (0:ht:tmax)'; Nt = length(t);

% Generate the initial condition for probability concentration, doing our
% best to replicate a delta function
c0 = zeros(length(x), 1);
zci = find(diff(sign(x)));  % find indices of zero crossings
if length(zci) == 1
    if abs(x(zci)) < abs(x(zci+1))
        c0(zci) = 1/hx;
    else
        c0(zci+1) = 1/hx;
    end
else
    c0(zci(2)) = 1/hx; 
end

% store probability concentrations for positive and negative drift
cp = crankNicolson_driftDiff(c0, x, t, p.D, 1);
cn = crankNicolson_driftDiff(c0, x, t, p.D, -1);

% store approximate FPT densities, one for each pair of threshold and drift
% direction
fMat = zeros(Nt, 4);
fMat(:,1) = cp(end-1,:)/hx;  % fHp
fMat(:,2) = cp(2,:)/hx;      % fLp
fMat(:,3) = cn(end-1,:)/hx;  % fHn
fMat(:,4) = cn(2,:)/hx;      % fLn

int_fMat = zeros(Nt,4); % from t to infinity
for n = Nt-1:-1:1
    for i = 1:4
        int_fMat(n,i) = int_fMat(n+1,i) + 1/2*(fMat(n,i)+fMat(n+1,i))*ht;
    end
end

% FPT densities conditioned on (WLOG) agent 1 deciding before agent 2
fHpcond = 2*fMat(:,1).*(int_fMat(:,1)+int_fMat(:,2));
fLpcond = 2*fMat(:,2).*(int_fMat(:,1)+int_fMat(:,2));
fHncond = 2*fMat(:,3).*(int_fMat(:,3)+int_fMat(:,4));
fLncond = 2*fMat(:,4).*(int_fMat(:,3)+int_fMat(:,4));

% Store the antiderivative of rho (see "intRho_x" function below). 
% Columns are 'crosses H given mu=1' and 'crosses -L given mu=-1'.
int_rhoMat = ones(Nt, 2); 

% pass Nt for dimension compatibility
int_rhoMat(:,1) = intRho_x(p, x, Nt, cp, 1);
int_rhoMat(:,2) = intRho_x(p, x, Nt, cn, -1);

%numSims = 5e4;
%figGenerator(p, intcpDenom, fMat(:,1:2), [fHpcond,fLpcond], t, tmax, numSims);

r = exp_reward(p, fHpcond, fLncond, int_rhoMat, t);
TF = exp_TF(fHpcond, fLpcond, fHncond, fLncond, t);
RR = 2*r / (TF + TI); 
end

function r = exp_reward(p, fHpcond, fLncond, int_rhoMat, t)
% Expected reward for each agent across all environments, WLOG agent 1.

probAtUpper = 1/2 * ( P_crossBefore(fHpcond, t) + ...
    P_instCrossAfter(fHpcond, int_rhoMat(:,1), t) );
probAtLower = 1/2 * ( P_crossBefore(fLncond, t) + ...
    P_instCrossAfter(fLncond, int_rhoMat(:,2), t) );

r = p.R1p/2 * probAtUpper + p.R1n/2 * probAtLower;
end

function TF = exp_TF(fHpcond, fLpcond, fHncond, fLncond, t)
% Expected time for the last agent to make a decision (WLOG for agent 2 to
% decide, conditioned on agent 1 deciding first).

TFp = trapz(t, t.*(fHpcond+fLpcond));
TFn = trapz(t, t.*(fHncond+fLncond));
TF = TFp/2 + TFn/2;
end

function P = P_crossBefore(fthreshcond, t)
% Probability that one agent crosses the decision threshold at "thresh",
% given that it decides before the other agent.

P = trapz(t, fthreshcond);
end

function P = P_instCrossAfter(fthreshcond, rhoVec, t)
% Probability that one agent is instantaneously "kicked" across the 
% decision threshold at "thresh", given the other agent's crossing of 
% "thresh" before.

fcondrho = fthreshcond.*rhoVec;
P = trapz(t, fcondrho);
end

function intRho = intRho_x(p, x, Nt, c, mu)
% Evaluate the x-antiderivative of rho (defined at the bottom of pg. 6 in
% my "Optimal Decision Strategies" notes) using the trapezoidal rule.

    intcnum = zeros(Nt, 1); intcdenom = intcnum;
    lowerBoundary = 1/2*(p.H1+p.L1)-p.qp1; 
    upperBoundary = 1/2*(p.H1+p.L1)+p.qn1;
    if sign(mu)==1 && lowerBoundary > p.L1
        idx = find(diff(sign( x - lowerBoundary )));
        if abs(x(idx(1)+1)-lowerBoundary) < abs(x(idx(1))-lowerBoundary)
            idx = idx + 1;
        end
        j = idx(1);
        for n = 1:Nt
            intcnum(n) = trapz(x(j:end), c(j:end,n));
            intcdenom(n) = trapz(x, c(:,n));
        end
    elseif sign(mu)==-1 && upperBoundary < p.H1
        idx = find(diff(sign( x - upperBoundary )));
        if abs(x(idx(1)+1)-upperBoundary) < abs(x(idx(1))-upperBoundary)
            idx = idx + 1;
        end
        j = idx(1);
        for n = 1:Nt
            intcnum(n) = trapz(x(1:j), c(1:j,n));
            intcdenom(n) = trapz(x, c(:,n));
        end
    else
        for n = 1:Nt
            intcnum(n) = trapz(x, c(:,n));
            intcdenom(n) = intcnum(n);
        end
    end
    
    intRho = intcnum ./ intcdenom;
end

function figGenerator(p, intcpDenom, fMatp, fMatpcond, t, tmax, numSims)
% create numerical vs Monte Carlo figures for
% (1) overall survival probability vs t;
% (2) unconditional FPT density (at H) vs t;
% (3) conditional FPT density (at H) vs t

% NOTE BELOW: using p assumes mu=1 for agent 1, mu=-1 for agent 2
w = tmax/200; % bin width
%p = parameters(p,'dt',min([10^(floor(log10(w))-1), 1e-3]));
[passTimesH, passTimesL] = passageTimes_oneAgent(p, numSims, 1);
[passTimes1H, passTimes1L, ~, ~] = passageTimes_twoAgent(p, numSims, 1, 1);

% (1) overall survival probability vs t
passTimes = sort([passTimesH passTimesL]);
survived = numSims*ones(1, int64(passTimes(end)/p.dt)+1);
for i = 1:length(passTimes)
    idx = int64(passTimes(i)/p.dt+1);
    survived(idx:end) = survived(idx:end) - 1;
end

figure;
plot(t, intcpDenom, '-r'); hold on
plot(0:p.dt:passTimes(end), survived/numSims, '.k', 'MarkerSize', 2);
title('Single-agent survival probability vs time, assuming \mu=1');
xlabel('Time'); ylabel('Survival probability')
legend('Crank-Nicolson', 'MC sims')


% (2) unconditional FPT density (at H) vs t;
bins = (w/2:w:tmax-w/2)';
FPTDH = zeros(length(bins),1);
FPTDL = zeros(length(bins),1);
for i=1:length(bins)
    FPTDH(i) = length(find(passTimesH >= w*(i-1) & passTimesH < w*i))/numSims/w;
    FPTDL(i) = length(find(passTimesL >= w*(i-1) & passTimesL < w*i))/numSims/w;
end

figure;
plot(t, fMatp(:,1), '-b'); hold on
plot(t, fMatp(:,2), '-r');
plot(bins, FPTDH, '.k', 'MarkerSize', 4);
plot(bins, FPTDL, '.k', 'MarkerSize', 4);
title('Single-agent FPT density at \pm threshold vs time, assuming \mu=1')
xlabel('Time')%,'FontSize',12)
ylabel('FPT Density')%,'FontSize',12)
legend('+ threshold', '- threshold', 'MC sims', '')


% (3) conditional FPT density (at H) vs t
numAgentOneFirst = length(passTimes1H)+length(passTimes1L);
FPTDHcond = zeros(length(bins),1);
FPTDLcond = zeros(length(bins),1);
for i=1:length(bins)
    FPTDHcond(i) = length(find(passTimes1H >= w*(i-1) & passTimes1H < w*i))/numAgentOneFirst/w;
    FPTDLcond(i) = length(find(passTimes1L >= w*(i-1) & passTimes1L < w*i))/numAgentOneFirst/w;
end

figure;
plot(t, fMatpcond(:,1), '-b'); hold on
plot(t, fMatpcond(:,2), '-r');
plot(bins, FPTDHcond, '.k', 'MarkerSize', 4);
plot(bins, FPTDLcond, '.k', 'MarkerSize', 4);
title('FPT density at \pm threshold (given agent 1 decides first) vs time, assuming \mu=1')
xlabel('Time')%,'FontSize',12)
ylabel('FPT Density')%,'FontSize',12)
legend('+ threshold', '- threshold', 'MC sims', '')

end