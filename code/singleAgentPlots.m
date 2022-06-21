% This script is a bit of a mess... need to come back and clean it up

global mu D
mu = 1; D = 1;

FPTfig = 0;
exitProbfig = 0;
RRfig = 0;
RRmaxfig = 1;
RR_expTfig = 0;
RR_exp_rfig = 0;

L = eps:.1:5;
H = L';
FPTmat = FPT(L, H);
condlExitMat = condlExitProb(L, H);
expT_RRMat = RR_expT(H,L);
exp_r_RRMat = RR_exp_r(H,L,1,1);

if FPTfig
    figure(1)
    h1 = contourf(L, H, FPTmat'); colorbar; hold on
    title('Mean Unconditional FPT (\mu=D=1)');
    xlabel('Lower threshold (magnitude)'); ylabel('Upper threshold');
    hold off
end

if exitProbfig
    figure(2)
    h2 = contourf(L, H, condlExitMat'); colorbar; hold on
    title('Probability of eventual exit through upper threshold H (\mu=D=1)');
    xlabel('Lower threshold (magnitude)'); ylabel('Upper threshold');
    hold off
end

if RRfig
    RH = 1; RL = 2; TI = 5;
    RRmat = zeros(length(H),length(L));
    %close all
    for i = 1:length(H)
        for j = 1:length(L)
            RRmat(i,j) = rewardRate([H(i) L(j)], RH, RL, TI);
        end
    end
    figure
    h3 = contourf(H, L, RRmat', 20);
    colorbar;
    title(['Reward Rate (R_L=' num2str(RL) ', R_H=' num2str(RH) ', T_I=' num2str(TI) ')']);
    xlabel('Upper threshold, H'); ylabel('Lower threshold, L');
    set(gcf,'color','w');
    hold off
end

if RRmaxfig
    TI = 5;
    Rp = 1:-.01:0; Rn = 1:.01:2;
    RRmat = zeros(length(Rp),1);
    X = zeros(length(Rp), 2);
    RRmatCN = zeros(length(Rp),1);
    XCN = zeros(length(Rp), 2);
    opts = optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-8, ...
        'StepTolerance',1e-8,'Algorithm','sqp','UseParallel',true,'DiffMaxChange',1);
    for i = 1:length(Rp)
        if i==1, X0 = [1.5 1.5]; 
        else, X0 = X(i-1,:); end
        problem = struct;
        RRfunc = @(X) -rewardRate(X, Rp(i), Rn(i), TI);
        problem.objective = RRfunc;
        problem.x0 = X0;
        problem.solver = 'fmincon';
        problem.lb = [eps eps];
        problem.ub = [10 10];
        problem.options = opts;
        %[params, RR] = fminsearchbnd( ...
        %    @(X) -RR_symmAgent_asymmThresh2(p, X, TI, N), ...%, LagRoots, LagWts), ...
        %    X0, [eps -inf eps eps], [inf -eps inf inf], opts);
        [params, RR] = fmincon(problem);
        RRmat(i) = -RR;
        X(i,:) = params;
        problem.objective = @(X) -RR_CrankNicolson(X, Rp(i), Rn(i), TI);
        [params, RR] = fmincon(problem);
        RRmatCN(i) = -RR;
        XCN(i,:) = params;
    end
    
    figure
    plot(Rp, RRmat, '-k'); hold on
    plot(Rp, RRmatCN, '--k');
    title('Optimized reward rate vs. R^+ (R^+ + R^- = 2)');
    xlabel('R^+'); ylabel('Optimal reward rate');
    set(gcf,'color','w');
    legend('RR', 'RR C-N')
    
    figure
    plot(Rp, X(:,1),'-b'); hold on
    plot(Rp, X(:,2),'-r');
    plot(Rp, XCN(:,1),'--b');
    plot(Rp, XCN(:,2),'--r');
    title('Optimized thresholds vs. R^+ (R^+ + R^- = 2)');
    xlabel('R^+'); ylabel('Optimal thresholds');
    set(gcf,'color','w');
    legend('+ threshold', '- threshold', '+ threshold C-N', '- threshold C-N')
end

if RR_expTfig
    figure
    contourf(L, H, expT_RRMat'); colorbar; hold on
    title('Expected Decision Time');
    xlabel('Lower threshold (magnitude)'); ylabel('Upper threshold');
    set(gcf,'color','w');
    hold off
end

if RR_exp_rfig
    figure
    contourf(L, H, exp_r_RRMat'); colorbar; hold on
    title('Expected Reward per Decision (R^+ = R^- = 1)');
    xlabel('Lower threshold (magnitude)'); ylabel('Upper threshold');
    set(gcf,'color','w');
    hold off
end

function expT = FPT(L,H)
global mu
    expT = (H.*phi(L+H) - (L+H).*phi(H) + L) ./ (phi(L+H)-1) ./ mu;
end

function pi_H = condlExitProb(L,H)
    pi_H = (phi(L)-1) ./ (phi(L+H)-1);
end

function val = phi(x)
global mu D
    val = exp(mu*x/D);
end

function RR = RR_CrankNicolson(X,RH,RL,TI)
global D

H = X(1); L = -X(2);
frac = H/(H-L);
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

hx = (H-L)/numInt;
min_char_t = 3*min([H, abs(L)]); tmax = 3*max([H, abs(L)]);
ht = min([min_char_t/1000, hx^2]); % adheres to CFL condition

x = (L:hx:H)';
t = (0:ht:tmax)'; Nt = length(t);

c0 = zeros(length(x), 1);
zci = find(diff(sign(x))); % find indices of zero crossings
if length(zci) == 1
    if abs(x(zci)) < abs(x(zci+1))
        c0(zci) = 1/hx;
    else
        c0(zci+1) = 1/hx;
    end
else
    c0(zci(2)) = 1/hx; 
end

cp = crankNicolson_driftDiff(c0, x, t, D, 1);
cn = crankNicolson_driftDiff(c0, x, t, D, -1);

fMat = zeros(Nt, 4);
fMat(:,1) = cp(end-1,:)/hx; % fHp
fMat(:,2) = cp(2,:)/hx; % fLp
fMat(:,3) = cn(end-1,:)/hx; % fHn
fMat(:,4) = cn(2,:)/hx; % fLn

r = RH/2 * trapz(t, fMat(:,1)) + RL/2 * trapz(t, fMat(:,4));
T = 1/2 * trapz(t, t.*(fMat(:,1)+fMat(:,2))) + ...
    1/2 * trapz(t, t.*(fMat(:,3)+fMat(:,4)));
RR = r / (T + TI); 
end

function RR = rewardRate(X,RH,RL,TI)
    H = X(1); L = X(2);
    RR = ( RH.*(1-exp(-L)) + RL.*(1-exp(-H)) ) ./ ...
        ( (1-exp(-L)).*(1-exp(-H)).*(L+H) + 2*TI.*(1-exp(-L-H)) );
end

function T = RR_expT(H,L)
    T = 1/2 * (1-exp(L)).*(1-exp(H)).*(L+H) ./ (exp(L+H)-1);
end

function r = RR_exp_r(H,L,Rp,Rn)
    r = 1/2 * ( Rp.*(exp(L)-1).*exp(H) + Rn.*(exp(H)-1).*exp(L) ) ...
        ./ (exp(L+H)-1);
end