function cfunc = c(p,x,t,agent,numImages,mu)
% Evaluates the Green's function expansion of an agent's probability
% concentration for the given inputs. See "Optimal Decision Strategies 
% for Agents with Asymmetric Decision Thresholds" (Fitzgerald) for more
% details.

    if agent == 1
        L = p.L1; H = p.H1;
    elseif agent == 2
        L = p.L2; H = p.H2;
    end
    cfunc = phi(p,x,t,mu);
    for n=1:numImages
        xi_L = xi(abs(L),H,n); xi_H = xi(H,abs(L),n);
        cfunc = cfunc + (-1)^n .* ( ...
            exp(-mu.*xi_L/p.D).*phi(p,x+2*xi_L,t,mu) + ...
            exp(mu.*xi_H/p.D).*phi(p,x-2*xi_H,t,mu) );
    end
    cfunc = cfunc ./ sqrt(4*pi*p.D.*t);
end

function fund_soln = phi(p,x,t,mu)
    fund_soln = exp(-(x-mu.*t).^2./(4*p.D.*t));
end

function val = xi(a,b,n)
    val = ceil(n/2)*a+floor(n/2)*b;
end