function M = Quadrupole(k,L)

% transfer matrix for a quadrupole (thick lens, first order approximation)
% taken from B. Holzer, Lattice Design in High-Energy Particle Accelerators
% CERN-2006-002, page 37

% we need just one case; the case with negative k will revert to the one
% with positive k

%sk = sign(k)*sqrt(abs(k));
if k~=0
    sk = sqrt(k);
    
    M = [cos(sk*L)      1/sk*sin(sk*L)  0              0                0 0 ;
        -sk*sin(sk*L)  cos(sk*L)       0              0                0 0;
        0              0               cosh(sk*L)     1/sk*sinh(sk*L)  0 0;
        0              0               sk*sinh(sk*L)  cosh(sk*L)       0 0;
        0              0                0              0               1 0;
        0              0                0              0               0 1];
else
    M = [1 L 0 0 0 0;
        0 1 0 0 0 0;
        0 0 1 L 0 0;
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1];
end