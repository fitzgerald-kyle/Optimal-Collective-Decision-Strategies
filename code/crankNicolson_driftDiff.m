function phiMat = crankNicolson_driftDiff(phi0, x, t, D, mu)
% assumes Dirichlet BCs
hx = x(2)-x(1); ht = t(2)-t(1);
Nx = length(x); Nt = length(t);
    
phiMat = zeros(Nx, Nt);
phiMat(:,1) = phi0;

for n = 2:Nt
    A = D*ht/(2*hx^2) + mu*ht/(4*hx); 
    B = 1 + D*ht/hx^2; 
    C = D*ht/(2*hx^2) - mu*ht/(4*hx); 
    E = zeros(Nx-1,1); F = E;
    for j = 2:Nx-1
        RHS = A*phiMat(j-1,n-1) + (2-B)*phiMat(j,n-1) + C*phiMat(j+1,n-1);
        E(j) = C/( B-A*E(j-1) );
        F(j) = ( RHS+A*F(j-1) )/( B-A*E(j-1) );
    end
    for j = Nx:-1:3
        phiMat(j-1,n) = E(j-1)*phiMat(j,n) + F(j-1);
    end
end
    
end