function [RAA,TCA] = ThreeSectionRT(Raa_1,Rbb_2,Rbb_1,Tba_1,Tcb_2,Tab_1,Tb)

% Reflection coefficients and phase change:
% Beam in logitudinal motion with 3 section discontinuity
% Breno Ebinuma Takiuti
% 28/08/2016

[dofb,~]=size(Rbb_1);

RAA = Raa_1 + Tab_1*pinv(eye(dofb)-Tb*Rbb_2*Tb*Rbb_1)*Tb*Rbb_2*Tb*Tba_1;
TCA = Tcb_2*pinv(eye(dofb)-Tb*Rbb_1*Tb*Rbb_2)*Tb*Tba_1;
