% Find the reflection and transmission from two sections modeled in FE
% using commercial package
% bend: both long and bending
% Breno Ebinuma Takiuti
% 27/04/2017

clear 
clc
close all


%% Beam Analytical solution

%% Material's Constants
% Material: Steel

rho=7800;     %mass per unit valume
E=2.06e11;  %Young's modulus

%% Geometric constants

b = 6e-3;              % base of the cross-section (tickness) (m)
ha = 18e-3;              % height of the cross-section (width) (m)
hb = 6e-3;
hc = ha;
La = 6e-4;
Lb = 6e-4;               % length of the element (x direction) (m)
Lc = La;
Sa = b*(ha);
Sb = b*hb;
Sc = Sa;
Ia = b*(ha)^3/12;
Ib = b*(hb)^3/12;
Ic = Ia;
beta_ab = Sb/Sa;
beta_bc = 1/beta_ab;

%% mass and stiffness from Ansys - plane
% K and M for the first section
% load KM_plane42_Lx005Ly0075_3el.mat
% load KM_plane42_Lx001Ly0075_12el.mat  %12 Plane elements with Lx = 0.001 Ly=0.0075
% load KM_plane43_Lx0006Ly009_15el.mat
% load KM_plane43_Lx0006Ly006_10el.mat
% load KM_plane43_Lx0006Ly012_20el.mat
load KM_plane43_Lx0006Ly018_30el.mat
% load KM_plane43_Lx0006Ly006_10el.mat
Ka0=full(K);
Ma0=full(M);

% Inod1 = 1:16;
% Inod2 = 9:24;
% Inod3 = 17:32;
% 
% Ka2 = zeros(2*size(Ka0));
% Ma2 = Ka2;
% 
% Ka2(Inod1,Inod1) = Ka0;
% Ka2(Inod2,Inod2) = Ka2(Inod2,Inod2)+Ka0;
% Ka2(Inod3,Inod3) = Ka2(Inod2,Inod2)+Ka0;
% Ma2(Inod1,Inod1) = Ma0;
% Ma2(Inod2,Inod2) = Ma2(Inod2,Inod2)+Ma0;
% Ma2(Inod3,Inod3) = Ma2(Inod2,Inod2)+Ma0;
% 
% Ka = Ka2;
% Ma = Ma2;
% 

Ka = Ka0;
Ma = Ma0;
Kc0 = Ka0;
Mc0 = Ma0;
Kc = Ka;
Mc = Ma;

% K and M for the second section
% load KM_plane42_Lx005Ly0025_1el.mat
load KM_plane43_Lx0006Ly006_10el.mat            % Cutoff 266kHz
% load KM_plane43_Lx0006Ly0072_12el.mat         % Cutoff 222kHz
% load KM_plane43_Lx0006Ly009_15el.mat          % Cutoff 170kHz
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
Kb=full(K);
Mb=full(M);

[ra, ca] = size(Ka0);
[rb, cb] = size(Kb);

Kabc = zeros(ra+ra,ca+ca);
Mabc = Kabc;

Kabc(1:ra,1:ca) = Ka0;
Mabc(1:ra,1:ca) = Ma0;

%% Boundary conditions when using the left eigenvectors
[ar,ac] = size(Ma);
[br,bc] = size(Mb);

Eb = zeros(ar/2,br/2);
Ca = zeros(ar/2);

% Describe the boundary (connection a-b)
I = eye(bc/2);
% Inoda = (ar/2-br/2)/2+1:(ar/2-br/2)/2+br/2; %SM
%     Inoda = (ar/2-44/2)/2+1:(ar/2-44/2)/2+br/2; %Quasi SM
%     Inoda = 1:br/2;                     %NS
%     Inoda = (ar/2-br/2)+1:ar/2;
Inoda = (ar-(ra/2)-br/2)/2+1:(ar-(ra/2)-br/2)/2+br/2; %SM2
Inodb = 1:br/2;

Inodab = [round(ra/2)+Inoda, ra+Inoda];
Inodbc = [ra+1:ra*2];

Ea = eye(ar/2);
Eb(Inoda,Inodb) = I;
Ec = Ea;
Ca(Inoda,Inoda) = I;
Cb = Eb;
Cc = Ca;

md = 1;

Kabc(Inodab,Inodab) = Kabc(Inodab,Inodab)+Kb;
Mabc(Inodab,Inodab) = Mabc(Inodab,Inodab)+Mb;
Kabc(Inodbc,Inodbc) = Kabc(Inodbc,Inodbc)+Kc0;
Mabc(Inodbc,Inodbc) = Mabc(Inodbc,Inodbc)+Mc0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2; %degree of freedom per node
na=30; %number of structural elements for the first section
nb=10; %number of elements for the second section
nc=na;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dofa = n*(na+1);       % number of degrees of freedom per side
dofb = n*(nb+1);
dofc = dofa;

%%
% Frequencies
fi = 100;
ff = 150000;
df = 1e3;
f = fi:df:ff;
% f = 10e3;

w = f(1)*2*pi;

%% Normalization
nor = 0;

%%

tol = 1e-4;
% Tried to change to descend sort for the complex numbers but did not work
% because the real and imaginary numbers are always descend, creating a
% discontinuity between complex and pure numbers.
[ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex6( w,Ka,Ma,La,nor,tol);
nmodes_a = length(kpai);
PsiQa_p0 = PsiQa_p(1:nmodes_a,:);
PsiQa_n0 = PsiQa_n(1:nmodes_a,:);
PsiFa_p0 = PsiFa_p(1:nmodes_a,:);
PsiFa_n0 = PsiFa_n(1:nmodes_a,:);

[ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex6( w,Kb,Mb,Lb,nor,tol);
nmodes_b = length(kpbi);
PsiQb_p0 = PsiQb_p(1:nmodes_b,:);
PsiQb_n0 = PsiQb_n(1:nmodes_b,:);
PsiFb_p0 = PsiFb_p(1:nmodes_b,:);
PsiFb_n0 = PsiFb_n(1:nmodes_b,:);

nmodes_c = nmodes_a;

w = 2*pi*f;
L = 1*Lb;
Lab = L+La;

for q=1:length(f)

    %% Waveguide A
    tol = 1e-2;
%     lim1 = 1;
%     lim2 = 1;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
      = PolySolve_complex6( w(q),Ka,Ma,La,nor,tol);
end

ai = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN

tolS = 1e-1;
% [a] = sortDiff(ai,f,tolS);
a = ai;

for q=1:length(f)
    %% B
%     tol = 1e-2;
%     lim1 = 0.1;
%     lim2 = 10;

%     [ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex( w,Kb,Mb,Lb,nor,tol);
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ] ...
      = PolySolve_complex6( w(q),Kb,Mb,Lb,nor,tol); %for 30-10 cases, use lim1 and lim2 default
end

bi = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
% [b] = sortDiff(bi,f);
b = bi;

for q=1:length(f)
    %% C
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
      = PolySolve_complex6( w(q),Kc,Mc,Lc,nor,tol);

end

ci = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
% [c] = sortDiff(ci,f);    
c = ci;

%% Should always use all propagating and nearfield waves, cant cut
%     mode_limita = 3;
    
    a.PhiQp = a.PhiQp(:,1:nmodes_a,:);
    a.PhiQn = a.PhiQn(:,1:nmodes_a,:);
    a.PsiQp = a.PsiQp(1:nmodes_a,:,:);
    a.PsiQn = a.PsiQn(1:nmodes_a,:,:);
    
    a.PhiFp = a.PhiFp(:,1:nmodes_a,:);
    a.PhiFn = a.PhiFn(:,1:nmodes_a,:);
    a.PsiFp = a.PsiFp(1:nmodes_a,:,:);
    a.PsiFn = a.PsiFn(1:nmodes_a,:,:);

    mode_limitc = nmodes_a;
    c.PhiQp = c.PhiQp(:,1:mode_limitc,:);
    c.PhiQn = c.PhiQn(:,1:mode_limitc,:);
    
    c.PsiQp = c.PsiQp(1:mode_limitc,:,:);
    c.PsiQn = c.PsiQn(1:mode_limitc,:,:);
    
    c.PhiFp = c.PhiFp(:,1:mode_limitc,:);
    c.PhiFn = c.PhiFn(:,1:mode_limitc,:);


for q=1:length(f)
    %% Analytical Solution
    
    [Maa,Kaa] = EB_Beam(rho,Sa,La,E,Ia);
    [Mab,Kab] = EB_Beam(rho,Sb,Lb,E,Ib);
    [Mac,Kac] = EB_Beam(rho,Sc,Lc,E,Ic);
    
    AA = eye(2);
    BB = AA;

    [RTa,PhiQ_pTa,PhiQ_nTa,PhiF_pTa,PhiF_nTa] = WM_reflection_beam(rho,Sa,E,Ia,w(q),AA,BB);
    [RTb,PhiQ_pTb,PhiQ_nTb,PhiF_pTb,PhiF_nTb] = WM_reflection_beam(rho,Sb,E,Ib,w(q),AA,BB);
    [RTc,PhiQ_pTc,PhiQ_nTc,PhiF_pTc,PhiF_nTc] = WM_reflection_beam(rho,Sc,E,Ic,w(q),AA,BB);
    
    kaa(q) = sqrt(w(q))*(rho*Sa/E/Ia)^(1/4);   % Wave number
    kbb(q) = sqrt(w(q))*(rho*Sb/E/Ib)^(1/4);   % Wave number
    kbc(q) = sqrt(w(q))*(rho*Sc/E/Ic)^(1/4); 
    
    lbpaa=exp(-1i*kaa*La);
    lbnaa=exp(1i*kaa*La);
    lbpab=exp(-1i*kbb*Lb);
    lbpac=exp(-1i*kbc*Lc);
    
    
    % Bending
    [RBTaa1(:,:,q),TBTba1(:,:,q)] = WA_reflection_beam_area(beta_ab,kaa(q),kbb(q));
    [RBTbb2(:,:,q),TBTcb2(:,:,q)] = WA_reflection_beam_area(beta_bc,kbb(q),kbc(q));
    
    RBTcc2(:,:,q) = RBTaa1(:,:,q);
    TBTbc2(:,:,q) = TBTba1(:,:,q);
    RBTbb1(:,:,q) = RBTbb2(:,:,q);
    TBTab1(:,:,q) = TBTcb2(:,:,q);
    
    % Bending transition matrix from interface 1 to 2
    TBTb = [exp(-1i*kbb(q)*L) 0; 0 exp(-kbb(q)*L)];
    
    % Bending scattering from 1 to 2
    [RBTAA(:,:,q),TBTCA(:,:,q)] = ThreeSectionRT(RBTaa1(:,:,q),RBTbb2(:,:,q),RBTbb1(:,:,q),TBTba1(:,:,q),TBTcb2(:,:,q),TBTab1(:,:,q),TBTb);

    % Longitudinal    
    kla(q) = (sqrt(rho/E)).*w(q); 
    klb(q) = kla(q); 
    klc(q) = kla(q); 
    
    [RLTaa1(q),TLTba1(q)] = WA_reflection_bar_area(beta_ab);
    [RLTbb2(q),TLTcb2(q)] = WA_reflection_bar_area(beta_bc);
    
    RLTcc2(q) = RLTaa1(q);
    TLTbc2(q) = TLTba1(q);
    RLTbb1(q) = RLTbb2(q);
    TLTab1(q) = TLTcb2(q);
    
    % Longitudinal transition matrix from interface 1 to 2
    TLTb = exp(-1i*klb(q)*L);
    
    % Longitudinal scattering from 1 to 2
    [RLTAA(:,:,q),TLTCA(:,:,q)] = ThreeSectionRT(RLTaa1(q),RLTbb2(q),RLTbb1(q),TLTba1(q),TLTcb2(q),TLTab1(q),TLTb);

    
    %% Numeric
    
    kPb = b.kp (1:nmodes_b,q);
    
    CP1ba = [PsiQa_p0*Ca*a.PhiQn(:,:,q) -PsiQa_p0*Cb*b.PhiQp(:,:,q);
        PsiFa_p0*Ea*a.PhiFn(:,:,q) -PsiFa_p0*Eb*b.PhiFp(:,:,q)];
    
    CP2ba = [-PsiQa_p0*Ca*a.PhiQp(:,:,q) PsiQa_p0*Cb*b.PhiQn(:,:,q);
        -PsiFa_p0*Ea*a.PhiFp(:,:,q) PsiFa_p0*Eb*b.PhiFn(:,:,q)]; 
  
    TRTba =(pinv(CP1ba)*CP2ba);
    [TRTr(q), TRTc(q)] = size(TRTba);
    TRT2Pba = TRTba;
    
    mode = round(TRTc(q)/2);
    
    RPaa1(:,:,q) = TRT2Pba(1:nmodes_a,1:nmodes_a);
    TPba1(:,:,q) = TRT2Pba(nmodes_a+1:end,1:nmodes_a);
    RPbb1(:,:,q) = TRT2Pba(nmodes_a+1:end,nmodes_a+1:end);
    TPab1(:,:,q) = TRT2Pba(1:nmodes_a,nmodes_a+1:end);
    
    % Bending transition matrix from interface 1 to 2
%     kPb = b.kp(:,q);
    TPb = diag(exp(-1i*kPb*L));
    
    %% Matrix inversion method as in Harland et al (2000) eq. (54) - b to c

    CP1cb = [PsiQa_p0*Cb*b.PhiQn(:,:,q) -PsiQa_p0*Cc*c.PhiQp(:,:,q);
        PsiFa_p0*Eb*b.PhiFn(:,:,q) -PsiFa_p0*Ec*c.PhiFp(:,:,q)];
    
    CP2cb = [-PsiQa_p0*Cb*b.PhiQp(:,:,q) PsiQa_p0*Cc*c.PhiQn(:,:,q);
        -PsiFa_p0*Eb*b.PhiFp(:,:,q) PsiFa_p0*Ec*c.PhiFn(:,:,q)];

    TRTPcb = (pinv(CP1cb)*CP2cb);
    
    RPbb2(:,:,q) = TRTPcb(1:nmodes_b,1:nmodes_b);
    TPcb2(:,:,q) = TRTPcb(nmodes_b+1:end,1:nmodes_b);
    RPcc2(:,:,q) = TRTPcb(nmodes_b+1:end,nmodes_b+1:end);
    TPbc2(:,:,q) = TRTPcb(1:nmodes_b,nmodes_b+1:end);
    
    % Scattering from a to c
    [RPAA(:,:,q),TPCA(:,:,q)] = ThreeSectionRT(RPaa1(:,:,q),RPbb2(:,:,q),RPbb1(:,:,q),TPba1(:,:,q),TPcb2(:,:,q),TPab1(:,:,q),TPb);

    %% Hybrid method by Renno
    
    Z1 = zeros(size(a.PhiQp(:,:,q)));
    Z2 = Z1';
    [rp,cp,~] = size(a.PhiQp);
    PhiQ_p = [a.PhiQp(:,:,q) Z1; Z1 c.PhiQp(:,:,q)];
    PhiQ_n = [a.PhiQn(:,:,q) Z1; Z1 c.PhiQn(:,:,q)];
    PhiF_p = [a.PhiFp(:,:,q) Z1; Z1 c.PhiFp(:,:,q)];
    PhiF_n = [a.PhiFn(:,:,q) Z1; Z1 c.PhiFn(:,:,q)];
    
    PsiQ_n = [a.PsiQn(:,:,q) Z2; Z2 c.PsiQn(:,:,q)];
    
    RA = eye(rp); RB = -eye(rp);
    Z3 = zeros(rp);
    
    R = [RA Z3; Z3 RB];
    
    InodE = [1:ra/2,ra*3/2+1:ra*2];
    InodI = ra/2+1:ra*3/2;
    
    D0 = (Kabc-w(q)^2*Mabc);
    
    %Method 1
%     DEE = D0(InodE,InodE);
%     DEI = D0(InodE,InodI);
%     DIE = D0(InodI,InodE);
%     DII = D0(InodI,InodI);

    %Method 2
    DTLL = D0(1:ra/2,1:ra/2);
    DTLR = D0(1:ra/2,ra*3/2+1:ra*2);
    DTRL = D0(ra*3/2+1:ra*2,1:rp);
    DTRR = D0(ra*3/2+1:ra*2,ra*3/2+1:ra*2);
    DTLO = D0(1:ra/2,ra/2+1:ra*3/2);
    DTRO = D0(ra*3/2+1:ra*2,ra/2+1:ra*3/2);
    DTOL = D0(ra/2+1:ra*3/2,1:ra/2);
    DTOR = D0(ra/2+1:ra*3/2,ra*3/2+1:ra*2);
    DTOO = D0(ra/2+1:ra*3/2,ra/2+1:ra*3/2);
    
    DLL = DTLL - DTLO*pinv(DTOO)*DTOL;
    DLR = DTLR - DTLO*pinv(DTOO)*DTOR;
    DRL = DTRL - DTRO*pinv(DTOO)*DTOL;
    DRR = DTRR - DTRO*pinv(DTOO)*DTOR;
    
    Dj = [DLL DLR; DRL DRR];    %Method 2
    
%     Dj = (DEE-DEI*DII^-1*DIE);    %Method 1
        
%     Srt(:,:,q) = pinv(PsiQ_n*(-Dj*R*PhiQ_n+R*PhiF_n))*PsiQ_n*(Dj*R*PhiQ_p-R*PhiF_p);
    Srt(:,:,q) = -pinv(Dj*R*PhiQ_n-R*PhiF_n)*(Dj*R*PhiQ_p-R*PhiF_p);

    RHAA(:,:,q) = Srt(1:3,1:3,q);
    THCA(:,:,q) = Srt(cp+1:cp+3,cp+1:cp+3,q);
    
    %%
%     
%      mode_limita = nmodes_a;
%     PhiQa_p2 = a.PhiQp(:,1:mode_limita,q);
%     PhiQa_n2 = a.PhiQn(:,1:mode_limita,q);
%     
%     PhiFa_p2 = a.PhiFp(:,1:mode_limita,q);
%     PhiFa_n2 = a.PhiFn(:,1:mode_limita,q);
%     
%     lpa2 = a.lp(1:mode_limita,q);
%     lna2 = a.ln(1:mode_limita,q);
% 
%     mode_limitb = nmodes_b;
%     PhiQb_p2 = b.PhiQp(:,1:mode_limitb,q);
%     PhiQb_n2 = b.PhiQn(:,1:mode_limitb,q);
%     
%     PsiQb_p2(:,:,q) = b.PsiQp(1:mode_limitb,:,q);
%     PsiQb_n2(:,:,q) = b.PsiQn(1:mode_limitb,:,q);
%     
%     PhiFb_p2 = b.PhiFp(:,1:mode_limitb,q);
%     PhiFb_n2 = b.PhiFn(:,1:mode_limitb,q);
%         
%     lpb2 = b.lp(1:mode_limitb,q);
%     lnb2 = b.ln(1:mode_limitb,q);
%     
%     mode_limitc = nmodes_c;
%     PhiQc_p2 = c.PhiQp(:,1:mode_limitc,q);
%     PhiQc_n2 = c.PhiQn(:,1:mode_limitc,q);
%     
%     PsiQc_p2(:,:,q) = c.PsiQp(1:mode_limitc,:,q);
%     PsiQc_n2(:,:,q) = c.PsiQn(1:mode_limitc,:,q);
%     
%     PhiFc_p2 = c.PhiFp(:,1:mode_limitc,q);
%     PhiFc_n2 = c.PhiFn(:,1:mode_limitc,q);
%         
%     lpc2 = c.lp(1:mode_limitc,q);
%     lnc2 = c.ln(1:mode_limitc,q);
%     
    
     %% Power matrix for the Plane element 1
    % MITROU (2015)
    Pa2 = (1i*w(q)/2)*[a.PhiQp(:,:,q)'*a.PhiFp(:,:,q) a.PhiQp(:,:,q)'*a.PhiFn(:,:,q);
        a.PhiQn(:,:,q)'*a.PhiFp(:,:,q) a.PhiQn(:,:,q)'*a.PhiFn(:,:,q)]-...
        [a.PhiFp(:,:,q)'*a.PhiQp(:,:,q) a.PhiFp(:,:,q)'*a.PhiQn(:,:,q);
        a.PhiFn(:,:,q)'*a.PhiQp(:,:,q) a.PhiFn(:,:,q)'*a.PhiQn(:,:,q)];
    Pc2 = (1i*w(q)/2)*[c.PhiQp(:,:,q)'*c.PhiFp(:,:,q) c.PhiQp(:,:,q)'*c.PhiFn(:,:,q);
        c.PhiQn(:,:,q)'*c.PhiFp(:,:,q) c.PhiQn(:,:,q)'*c.PhiFn(:,:,q)]-...
        [c.PhiFp(:,:,q)'*c.PhiQp(:,:,q) c.PhiFp(:,:,q)'*c.PhiQn(:,:,q);
        c.PhiFn(:,:,q)'*c.PhiQp(:,:,q) c.PhiFn(:,:,q)'*c.PhiQn(:,:,q)];
    

    PrPP2(q) = abs(RPAA(1,1,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(1,1));
    PrPL2(q) = abs(RPAA(2,1,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(1,1));
    PrPN2(q) = abs(RPAA(3,1,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(1,1));
    PtPP2(q) = abs(TPCA(1,1,q))^2*(Pc2(1,1)/Pa2(1,1));
    PtPL2(q) = abs(TPCA(2,1,q))^2*(Pc2(2,2)/Pa2(1,1));
    PtPN2(q) = abs(TPCA(3,1,q))^2*(Pc2(3,3)/Pa2(1,1));
  
    PrLP2(q) = abs(RPAA(1,2,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(2,2));
    PrLL2(q) = abs(RPAA(2,2,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(2,2));
    PrLN2(q) = abs(RPAA(3,2,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(2,2));

    PtLP2(q) = abs(TPCA(1,2,q))^2*(Pc2(1,1)/Pa2(2,2));
    PtLL2(q) = abs(TPCA(2,2,q))^2*(Pc2(2,2)/Pa2(2,2));
    PtLN2(q) = abs(TPCA(3,2,q))^2*(Pc2(3,3)/Pa2(2,2));

    
    %% Kinetic Energy - BENDING INCIDENT
    
    AMP_I = zeros(nmodes_a,1);
    AMP_I(1) = 1;
    AMP_R = RPAA(:,:,q)*AMP_I;
    AMP_T = TPCA(:,:,q)*AMP_I;
    
    % Bending Propagating Positive in A  (incident)
    [EIbbkin_avL(q),EIbbpot_av(q),EIbbpot_perx(q),EIbbpot_pery(q),PIbbp(q),PRbbx(q),PIbby(q),EIbbkin_avLx(q),EIbbkin_avLy(q),veIbbp(q)] = KineticEn3(a.PhiQp(:,:,q),a.lp(q),Ma,Ka,w(q),La,(AMP_I));
    
    % Bending Propagating Negative in A  (reflected)
    [ERbbkin_avL(q),ERbbpot_av(q),ERbbpot_perx(q),ERbbpot_pery(q),PRbbn(q),PRbbx(q),PRbby(q),ERbbkin_avLx(q),ERbbkin_avLy(q),veRbbn(q)] = KineticEn3(a.PhiQn(:,:,q),a.ln(q),Ma,Ka,w(q),La,(AMP_R));

    % Bending Propagating Positive in B (transmitted) 
    [ETbbkin_avL(q),ETbbpot_av(q),ETbbpot_perx(q),ETbbpot_pery(q),PTbbp(q),PTbbx(q),PTbby(q),ETbbkin_avLx(q),ETbbkin_avLy(q),veTbbp(q)]  = KineticEn3(c.PhiQp(:,:,q),c.lp(q),Mc,Kc,w(q),Lc,(AMP_T));
   
    %% Kinetic Energy - BENDING INCIDENT Analytical
    
    [nmodes_aa,~] = size(RBTAA);
    AMP_Ia = zeros(nmodes_aa,1);
    AMP_Ia(1) = 1;
    AMP_Ra = RBTAA(:,:,q)*AMP_Ia;
    AMP_Ta = TBTCA(:,:,q)*AMP_Ia;
    
    % Bending Propagating Positive in A  (incident)
    [EIbbkin_avLa(q),~,~,~,~,~,~,EIbbkin_avLxa(q),EIbbkin_avLya(q),~] = KineticEn3(PhiQ_pTa,lbpaa(q),Maa,Kaa,w(q),La,(AMP_Ia));
    
    % Bending Propagating Negative in A  (reflected)
    [ERbbkin_avLa(q),~,~,~,~,~,~,ERbbkin_avLxa(q),ERbbkin_avLya(q),~] = KineticEn3(PhiQ_nTa,lbnaa(q),Maa,Kaa,w(q),La,(AMP_Ra));

    % Bending Propagating Positive in B (transmitted) 
    [ETbbkin_avLa(q),~,~,~,~,~,~,ETbbkin_avLxa(q),ETbbkin_avLya(q),~]  = KineticEn3(PhiQ_pTc,lbpac(q),Mac,Kac,w(q),Lc,(AMP_Ta));

    %% Kinetic Energy - LONGITUDINAL INCIDENT 
 
    AMP_I = zeros(nmodes_a,1);
    AMP_I(2) = 1;
    AMP_R = RPAA(:,:,q)*AMP_I;
    AMP_T = TPCA(:,:,q)*AMP_I;
    
    % Bending Propagating Positive in A  (incident)
    [EIllkin_avL(q),EIllpot_av(q),EIllpot_perx(q),EIllpot_pery(q),PIlbp(q),PRllx(q),PIlly(q),EIllkin_avLx(q),EIllkin_avLy(q),veIllp(q)] = KineticEn3(a.PhiQp(:,:,q),a.lp(q),Ma,Ka,w(q),La,AMP_I);
    
    % Bending Propagating Negative in A  (reflected)
    [ERlbkin_avL(q),ERlbpot_av(q),ERlbpot_perx(q),ERlbpot_pery(q),PRbbn(q),PRlbx(q),PRlby(q),ERlbkin_avLx(q),ERlbkin_avLy(q),veRlbn(q)] = KineticEn3(a.PhiQn(:,:,q),a.ln(q),Ma,Ka,w(q),La,AMP_R);

    % Bending Propagating Positive in C (transmitted) 
    [ETlbkin_avL(q),ETlbpot_av(q),ETlbpot_perx(q),ETlbpot_pery(q),PTbbp(q),PTlbx(q),PTlby(q),ETlbkin_avLx(q),ETlbkin_avLy(q),veTlbp(q)]  = KineticEn3(c.PhiQp(:,:,q),c.lp(q),Mc,Kc,w(q),Lc,AMP_T);

end

TRTP = [RPAA; TPCA];

% Analyticlal Reflection and Transmission
R_WM(1,:) = RBTAA(1,1,:);
T_WM(1,:) = TBTCA(1,1,:);
R_WM(2,:) = RLTAA;
T_WM(2,:) = TLTCA;
R_WM(3,:) = RBTAA(2,1,:);
T_WM(3,:) = TBTCA(2,1,:);

% PLANE1 Extract the R&T from scaterring matrix
R_WFE2(1,:) = reshape(RPAA(1,1,:),[1 length(f)]);
R_WFE2(2,:) = reshape(RPAA(2,2,:),[1 length(f)]);
T_WFE2(1,:) = reshape(TPCA(1,1,:),[1 length(f)]);
T_WFE2(2,:) = reshape(TPCA(2,2,:),[1 length(f)]);

R_WFE1(1,:) = reshape(RPaa1(1,1,:),[1 length(f)]);
R_WFE1(2,:) = reshape(RPaa1(2,2,:),[1 length(f)]);
T_WFE1(1,:) = reshape(TPab1(1,1,:),[1 length(f)]);
T_WFE1(2,:) = reshape(TPba1(2,2,:),[1 length(f)]);

% Cutoff frequencies
Cutoff_ai = [88e3 154e3 166e3 177e3 266e3]; 
Cutoff_bi = [132e3 231e3 249e3 268e3];                           % b12
% Cutoff_bi = 266e3;

Cutoff_a = round(Cutoff_ai/df+1); 
Cutoff_b = round(Cutoff_bi/df);
% Cutoff_a = round(([156e3 166e3 ]-fi)/df); 
% Cutoff_b = round((266e3-fi)/df);


 %% Cutoff plots
% Use the plotCutOffs only for the non zero coefficients, the others plot
% them normally
% 
figure()
plot(f,abs(reshape(RPAA(1,1,:),[1 length(f)])),'b','LineWidth',2);
hold on
plot(f,abs(reshape(RPAA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(RPAA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(RPAA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPAA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);

plot(f,abs(reshape(RBTAA(1,1,:),[1 length(f)])),'b:')
plot(f,abs(reshape(RBTAA(2,1,:),[1 length(f)])),'b:')

plot(f,abs(reshape(RHAA(1,1,:),[1 length(f)])),'r--')
plot(f,abs(reshape(RHAA(2,1,:),[1 length(f)])),'r--')

for ii = 1:length(Cutoff_a)
    plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
end
for ii = 1:length(Cutoff_b)
    plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
end
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
 axis([fi ff 0 2])

figure()
plot(f,abs(reshape(TPCA(1,1,:),[1 length(f)])),'b','LineWidth',2);
hold on
plot(f,abs(reshape(TPCA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(TPCA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(TPCA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPCA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);

plot(f,abs(reshape(TBTCA(1,1,:),[1 length(f)])),'b:')
plot(f,abs(reshape(TBTCA(2,1,:),[1 length(f)])),'b:')

plot(f,abs(reshape(THCA(1,1,:),[1 length(f)])),'r--')
plot(f,abs(reshape(THCA(2,1,:),[1 length(f)])),'r--')

plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
 axis([fi ff 0 2])
 
%  figure()
% plot(f,abs(reshape(RPAA(1,2,:),[1 length(f)])),'b','LineWidth',2);
% hold on
% plot(f,abs(reshape(RPAA(2,2,:),[1 length(f)])),'b','LineWidth',2);
% % plotCutOffs(f,abs(reshape(RPAA(3,2,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% % plotCutOffs(f,abs(reshape(RPAA(4,2,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% % plotCutOffs(f,abs(reshape(RPAA(5,2,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);
% 
% plot(f,abs(R_WM(2,:,:)),'b:')
% 
% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])
% 
% figure()
% plot(f,abs(reshape(TPCA(1,2,:),[1 length(f)])),'b','LineWidth',2);
% hold on
% plot(f,abs(reshape(TPCA(2,2,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(TPCA(3,2,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(TPCA(4,2,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPCA(5,2,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);
% 
% plot(f,abs(T_WM(2,:)),'b:')
% 
% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])
%  
 
%% Energy plots
% 
%  % Bending Incident TA Kinetic Energy
%  figure()
%  plot(f,abs(EIbbkin_avLy./EIbbkin_avL)*100,'b*','LineWidth',1)
%  hold on
%  plot(f,abs(EIbbkin_avLx./EIbbkin_avL)*100,'rp','LineWidth',1)
%  plot(f,abs(EIbbkin_avLya./EIbbkin_avLa)*100,'k*','LineWidth',1)
%  plot(f,abs(EIbbkin_avLxa./EIbbkin_avLa)*100,'kp','LineWidth',1)
% 
%  legend('kEy_{inc}/kE_{inc}','kEx_{inc}/kE_{inc}')
% % %  title('Symmetric discontinuity')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Time-averaged kinetic energy [%]','FontName','Times New Roman','FontSize',12)
%  set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi f(q) -1 101])
% 
% 
% for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
%  figure()
%  plot(f,(ERbbkin_avLy./ERbbkin_avL)*100,'b*','LineWidth',1)
%  hold on
%  plot(f,(ERbbkin_avLx./ERbbkin_avL)*100,'rp','LineWidth',1)
%  plot(f,(ETbbkin_avLy./ETbbkin_avL)*100,'go','LineWidth',1)
%  plot(f,(ETbbkin_avLx./ETbbkin_avL)*100,'ms','LineWidth',1)
% 
%  legend('kEy_{refl}/kE_{refl}','kEx_{refl}/kE_{refl}','kEy_{trans}/kE_{trans}','kEx_{trans}/kE_{trans}')
% % %  title('Symmetric discontinuity')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Time-averaged kinetic energy [%]','FontName','Times New Roman','FontSize',12)
%  set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi f(q) -1 101])
% for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
