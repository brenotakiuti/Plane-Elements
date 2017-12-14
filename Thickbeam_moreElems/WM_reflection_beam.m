function [RT,PhiQ_pT,PhiQ_nT,PhiF_pT,PhiF_nT] = WM_reflection_beam(rho,A,E,In,w,AA,BB)

k = abs(sqrt(w)*(rho*A/E/In)^(1/4));

PhiQ_pT = [1 1; -1i*k -k];

PhiQ_nT = [1 1; 1i*k k];

PhiF_pT = E*In*k^2*[1i*k -k; -1 1];

PhiF_nT = E*In*k^2*[-1i*k k; -1 1];
size(PhiQ_pT);
% Use the analytical matrices
Ct = AA*PhiQ_nT - BB*PhiF_nT;
Dt = BB*PhiF_pT - AA*PhiQ_pT;

RT = (inv(Ct)*Dt);
