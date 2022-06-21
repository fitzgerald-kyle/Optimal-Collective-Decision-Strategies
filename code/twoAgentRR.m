% Two-agent reward rate plots

% Symmetric agent, asymmetric thresholds. 
TI = 5;
N = 1e4; % image terms

% See each "if" block for figure descriptions.
fig1 = 1;
fig2 = 0;
fig3 = 0;

if fig1
% Fix rewards and kicks, sweep over thresholds. Plot RR.
    H = linspace(1.6, 2, 15); L = -linspace(1.6, 2, 15); Rp = 1; Rn = 1;
    RRmat = zeros(length(H),length(L));
    for i= 1:length(H)
        for j= 1:length(L)
            qp = inf;%H(i); 
            qn = inf;
            p = parameters('qp1',qp,'qn1',qn,'H1',H(i),'L1',L(j),'R1p',Rp,'R1n',Rn);
            %tic;
            RRmat(i,j) = RR_symmAgent_asymmThresh_simple3(p, [H(i) L(j) qp qn], TI);
            %T=toc;
        end
        disp(i/length(H)*100);
    end
    
    figure
    contourf(H, L, RRmat', 20); colorbar;
    title(['Reward Rate (R^-=' num2str(Rn) ', R^+=' num2str(Rp) ...
        ', T_I=' num2str(TI) ')']);
    xlabel('+ threshold, H'); ylabel('- threshold, L');
    set(gcf,'color','w');
    hold off
end

if fig2
% Fix rewards and thresholds, sweep over kicks. Plot RR.

    qp = linspace(.1, 1.5, 15); qn = qp; Rp = .4; Rn = 1;
    RRmat = zeros(length(qp),length(qn));
    for i= 1:length(qp)
        for j= 1:length(qn)
            H = 1.05; 
            L = -.35;
            p = parameters('qp1',qp(i),'qn1',qn(j),'H1',H,'L1',L,'R1p',Rp,'R1n',Rn);
            %tic;
            RRmat(i,j) = RR_symmAgent_asymmThresh(p, [H L qp(i) qn(j)], TI, N);
            %T=toc;
        end
    end
    
    figure
    contourf(qp, qn, RRmat', 20); colorbar;
    title(['Reward Rate (R^-=' num2str(Rn) ', R^+=' num2str(Rp) ...
        ', T_I=' num2str(TI) ')']);
    xlabel('+ kick'); ylabel('- kick');
    set(gcf,'color','w');
    hold off
end

if fig3
% Maximize RR wrt thresholds and kicks while sweeping over rewards,
% adhering to Rp + Rn = 2. Uses fmincon.

    Rp = 1:-.01:0;%linspace(.05,1,20); 
    Rn = 1:.01:2;
    RRmat = zeros(1,length(Rp));
    RRsim_mat = zeros(1,length(Rp));
    RRsim_std = zeros(1,length(Rp));
    %RRsim_Q3 = zeros(1,length(Rp));
    X = zeros(length(Rp),4);%4);
    opts = optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-8, ...
        'StepTolerance',1e-8,'FiniteDifferenceType','central','Algorithm','sqp',...
        'UseParallel',true,'DiffMaxChange',1);
    %opts = optimset('Display','iter','TolFun',1e-8, ...
    %    'TolX',1e-8,'FinDiffType','central','MaxIter',100);
    LB = 1e-2; UB = 10;
    for i= 1:length(Rp)
        if i==1, X0 = [1.68 -1.68 1.68 1.68]; 
        else, X0 = X(i-1,:); end
        p = parameters('R1p',Rp(i),'R1n',Rn(i),'R2p',Rp(i),'R2n',Rn(i));
        %warning('off')
        problem = struct;
        problem.objective = @(X) -RR_symmAgent_asymmThresh_simple3(p, X, TI);%, N);
        problem.x0 = X0;
        problem.solver = 'fmincon';
        problem.lb = [LB -UB 0 0];
        problem.ub = [UB -LB 2*(LB+UB) 2*(LB+UB)];
        problem.options = opts;
        problem.Aineq = [-1,1,1,0; -1,1,0,1];%[x(3) - (x(1) - x(2)), x(4) - (x(1) - x(2))]
        problem.bineq = [0;0];
        %problem.nonlcon = @kickcon;
        [params, RR] = fmincon(problem);
        RRmat(i) = -RR;
        H = params(1); L = params(2); qp = params(3); qn = params(4);
        [RRsim_mat(i), RRsim_std(i)] = ...
            RRsim_simple(parameters('p', p, 'H1', H, 'L1', L, 'qp1', qp, 'qn1', qn, ...
            'H2', H, 'L2', L, 'qp2', qp, 'qn2', qn), TI);
        X(i,:) = params;
    end
    
    figure
    plot(Rp, RRmat, '-b'); hold on
    plot(Rp, RRsim_mat, '-r')
    p1 = fill([Rp fliplr(Rp)], ...
        [RRsim_mat-RRsim_std fliplr(RRsim_mat+RRsim_std)], 'r', 'FaceAlpha', 0.2);
    title('Optimized reward rate vs. positive reward (R^+ + R^- = 2)');
    xlabel('Positive reward, R^+'); ylabel('Optimal reward rate');
    legend('theoretical', 'sims')
    set(gcf,'color','w');
    
    figure
    plot(Rp, X(:,1),'-b'); hold on
    plot(Rp, X(:,2),'-r');
    title('Optimized thresholds vs. positive reward (R^+ + R^- = 2)');
    xlabel('Positive reward, R^+'); ylabel('Optimal thresholds');
    legend('H', '-L')
    set(gcf,'color','w');
    
    figure
    plot(Rp, X(:,3),'-b'); hold on
    plot(Rp, X(:,4),'-r');
    title('Optimized kicks vs. positive reward (R^+ + R^- = 2)');
    xlabel('Positive reward, R^+'); ylabel('Optimal kicks');
    legend('q_+', 'q_-')
    set(gcf,'color','w');
    hold off
end

%function [c,ceq] = kickcon(x)
%    c = [x(3) - (x(1) - x(2)), x(4) - (x(1) - x(2))];
%    ceq = [];
%end

function [RR,RRstd] = RRsim_simple(p, TI)
% Monte Carlo simulated reward rates

    n = 5e3;
    tvec = zeros(1,n);
    rvec = zeros(1,n);
    for i = 1:n
        mu = sign(2*rand-1);
        [X1,~,over1,over2,~] = ...
            diffusionTrial_twoAgent_simple(parameters(p,'mu1',mu,'mu2',mu));
        tvec(i) = (length(X1)-1)*p.dt;
        if mu == 1
            rvec(i) = p.R1p*over1 + p.R2p*over2;
        elseif mu == -1
            rvec(i) = p.R1n*(~over1) + p.R2n*(~over2);
        end
    end
    rbar = mean(rvec); tbar = mean(tvec);
    RR = rbar / (tbar+TI);
    %RRQ1 = quantile(rvec./(tvec+TI), 0.25);
    %RRQ3 = quantile(rvec./(tvec+TI), 0.75);
    rstd = std(rvec); tstd = std(tvec);
    RRstd = sqrt(rstd^2/(tbar+TI)^2 + (rbar*tstd)^2/(tbar+TI)^4);
end