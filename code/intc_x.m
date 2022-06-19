function intcfunc = intc_x(p,x,t,agent,numImages,mu)
% Evaluates the x-antiderivative of the Green's function expansion of an 
% agent's probability concentration for the given inputs. See "Optimal 
% Decision Strategies for Agents with Asymmetric Decision Thresholds" 
% (Fitzgerald) for more details.

    if agent == 1
        L = p.L1; H = p.H1;
    elseif agent == 2
        L = p.L2; H = p.H2;
    end
    intcfunc = intphi(p,x,t,mu);
    for n=1:numImages
        xi_L = xi(abs(L),H,n); xi_H = xi(H,abs(L),n);
        intcfunc = intcfunc + (-1)^n .* ( ...
            exp(-mu.*xi_L./p.D).*intphi(p,x+2*xi_L,t,mu) + ...
            exp(mu.*xi_H./p.D).*intphi(p,x-2*xi_H,t,mu) );
    end
end

function f = intphi(p,x,t,mu)
    f = 1/2 * erf( (x-mu.*t)./(2*sqrt(p.D.*t)) );
end

function val = xi(a,b,n)
    val = ceil(n/2)*a+floor(n/2)*b;
end