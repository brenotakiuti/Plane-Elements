% FAIL
function  [RPAA,TPCA] = AreaChangePlaneWFE(w,beta)
% Function for calculating the three section solution as function of the
% non dimensional thickness
% Breno Ebinuma Takiuti
% 04/07/2017

%% Geometric constants
wa = 6e-3;              % base of the cross-section (tickness) (m)
Lb = 6e-4;

%% mass and stiffness from Ansys - plane
% K and M for the first section
load KM_plane42_Lx001Ly0075_12el.mat  %12 Plane elements with Lx = 0.001 Ly=0.0075
% load KM_plane43_Lx0006Ly006_10el.mat
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
Ka=full(K);
Ma=full(M);

Kc = Ka;
Mc = Ma;

% K and M for the second section
load KM_plane42_Lx001Ly0025_4el.mat  %4 Plane elements with Lx=0.001 Ly=0.0025
% load KM_plane43_Lx0006Ly006_10el.mat
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
Kb=full(K);
Mb=full(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2; %degree of freedom per node
na=30; %number of structural elements for the first section
nb=10; %number of elements for the second section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dofa = n*(na+1);       % number of degrees of freedom per side

%%
% Frequencies
wi = 1e3*2*pi;

%% Normalization
nor = 0;

%%

tol = 1e-3;

[ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex( wi,Ka,Ma,Le,nor,tol);
nmodes_a = length(knai);
% mode_limita = 5;
PsiQa_p0 = PsiQa_p(1:nmodes_a,:);
PsiQa_n0 = PsiQa_n(1:nmodes_a,:);
PsiFa_p0 = PsiFa_p(1:nmodes_a,:);
PsiFa_n0 = PsiFa_n(1:nmodes_a,:);

[ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex( wi,Kb,Mb,Lb,nor,tol);
nmodes_b = length(knbi);
PsiQb_p0 = PsiQb_p(1:nmodes_b,:);
PsiQb_n0 = PsiQb_n(1:nmodes_b,:);
PsiFb_p0 = PsiFb_p(1:nmodes_b,:);
PsiFb_n0 = PsiFb_n(1:nmodes_b,:);

nmodes_c = nmodes_a;


%% Waveguide A

tol = 1e-4;

[ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex( w,Ka,Ma,Le,nor,tol,3);

%% B
[ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex( w,Kb,Mb,Lb,nor,tol,3);

%% C
[ PhiQc_p1,PhiQc_n1,PhiFc_p1,PhiFc_n1,PsiQc_p1,PsiFc_p1,PsiQc_n1,PsiFc_n1,lpc1,lnc1,sc,kpci,knci ] = PolySolve_complex( w,Kc,Mc,Lc,nor,tol,3);

%% Should always use all propagating and nearfield waves, cant cut
%     mode_limita = 7;

PhiQa_p1 = PhiQa_p(:,1:nmodes_a);
PhiQa_n1 = PhiQa_n(:,1:nmodes_a);

PhiFa_p1 = PhiFa_p(:,1:nmodes_a);
PhiFa_n1 = PhiFa_n(:,1:nmodes_a);

%     mode_limitb = 5;
mode_limitb = length(knbi);
PhiQb_p1 = PhiQb_p(:,1:mode_limitb);
PhiQb_n1 = PhiQb_n(:,1:mode_limitb);

PsiQb_p1(:,:,q) = PsiQb_p(1:mode_limitb,:);
PsiQb_n1(:,:,q) = PsiQb_n(1:mode_limitb,:);

PhiFb_p1 = PhiFb_p(:,1:mode_limitb);
PhiFb_n1 = PhiFb_n(:,1:mode_limitb);

%     PsiFb_p1(:,:,q) = PsiFb_p(1:mode_limitb,:);
%     PsiFb_n1(:,:,q) = PsiFb_n(1:mode_limitb,:);

%% Boundary conditions when using the left eigenvectors
[ar,ac] = size(Ma);
[br,bc] = size(Mb);

Eb = zeros(ar/2,br/2);
Ca = zeros(ar/2);

% Describe the boundary (connection a-b)
I = eye(bc/2);
%     Inoda = 21:42; %SM
Inoda = 9:18; %SM
%     Inoda = 1:br/2; %NS
Inodb = 1:br/2;

Ea = eye(ar/2);
Eb(Inoda,Inodb) = I;
Ec = Ea;
Ca(Inoda,Inoda) = I;
Cb = Eb;
Cc = Ca;

md = 1;

CP1ba = [PsiQa_p0*Ca*PhiQa_n1 -PsiQa_p0*Cb*PhiQb_p1;
    PsiFa_p0*Ea*PhiFa_n1 -PsiFa_p0*Eb*PhiFb_p1];

CP2ba = [-PsiQa_p0*Ca*PhiQa_p1 PsiQa_p0*Cb*PhiQb_n1;
    -PsiFa_p0*Ea*PhiFa_p1 PsiFa_p0*Eb*PhiFb_n1];

TRTba =(pinv(CP1ba)*CP2ba);
[TRTr(q), TRTc(q)] = size(TRTba);
TRT2Pba(:,:,q) = TRTba;

mode = round(TRTc(q)/2);

RPaa1 = TRT2Pba(1:nmodes_a,1:nmodes_a,q);
TPba1 = TRT2Pba(nmodes_a+1:end,1:nmodes_a,q);
RPbb1 = TRT2Pba(nmodes_a+1:end,nmodes_a+1:end,q);
TPab1 = TRT2Pba(1:nmodes_a,nmodes_a+1:end,q);

% Bending transition matrix from interface 1 to 2
TPb = diag(exp(-1i*kPb*Lb));

% Matrix inversion method as in Harland et al (2000) eq. (54) - b to c
CP1cb = [PsiQa_p0*Cb*PhiQb_n1 -PsiQa_p0*Cc*PhiQc_p1;
    PsiFa_p0*Eb*PhiFb_n1 -PsiFa_p0*Ec*PhiFc_p1];

CP2cb = [-PsiQa_p0*Cb*PhiQb_p1 PsiQa_p0*Cc*PhiQc_n1;
    -PsiFa_p0*Eb*PhiFb_p1 PsiFa_p0*Ec*PhiFc_n1];

TRTPcb(:,:,q) = (pinv(CP1cb)*CP2cb);

RPbb2 = TRTPcb(1:nmodes_b,1:nmodes_b,q);
TPcb2 = TRTPcb(nmodes_b+1:end,1:nmodes_b,q);
RPcc2 = TRTPcb(nmodes_b+1:end,nmodes_b+1:end,q);
TPbc2 = TRTPcb(1:nmodes_b,nmodes_b+1:end,q);

% Scattering from a to c
[RPAA,TPCA] = ThreeSectionRT(RPaa1,RPbb2,RPbb1,TPba1,TPcb2,TPab1,TPb);

