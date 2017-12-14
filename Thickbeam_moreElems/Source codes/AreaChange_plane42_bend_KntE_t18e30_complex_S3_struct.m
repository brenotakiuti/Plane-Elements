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
load KM_plane42_Lx001Ly0075_12el.mat  %12 Plane elements with Lx = 0.001 Ly=0.0075
% load KM_plane43_Lx0006Ly009_15el.mat
% load KM_plane43_Lx0006Ly006_10el.mat
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
% load KM_plane43_Lx0006Ly006_10el.mat
Ka=full(K);
Ma=full(M);

Kc = Ka;
Mc = Ma;

% K and M for the second section
load KM_plane42_Lx001Ly0025_4el.mat
% load KM_plane43_Lx0006Ly006_10el.mat            % Cutoff 266kHz
% load KM_plane43_Lx0006Ly0072_12el.mat         % Cutoff 222kHz
% load KM_plane43_Lx0006Ly009_15el.mat          % Cutoff 170kHz
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
Kb=full(K);
Mb=full(M);

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
ff = 250000;
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
[ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex4( w,Ka,Ma,La,nor,tol);
nmodes_a = length(kpai);
% nmodes_a = 3;
PsiQa_p0 = PsiQa_p(1:nmodes_a,:);
PsiQa_n0 = PsiQa_n(1:nmodes_a,:);
PsiFa_p0 = PsiFa_p(1:nmodes_a,:);
PsiFa_n0 = PsiFa_n(1:nmodes_a,:);

[ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex4( w,Kb,Mb,Lb,nor,tol);
nmodes_b = length(kpbi);
nmodes_b = 3;
% nmodes_b = 2;
PsiQb_p0 = PsiQb_p(1:nmodes_b,:);
PsiQb_n0 = PsiQb_n(1:nmodes_b,:);
PsiFb_p0 = PsiFb_p(1:nmodes_b,:);
PsiFb_n0 = PsiFb_n(1:nmodes_b,:);

nmodes_c = nmodes_a;

w = 2*pi*f;
% L = Lb;
% L = 0.18;
L = 0.09;
% a = struct('PhiQp','PhiQn','PhiFp','PhiFn','PsiQp','PsiFp','PsiQn','PsiFn','lp','ln','s','kp','kn','PureN');

for q=1:length(f)

    %% Waveguide A
    %30-10
%     tol = 1e-4;
%     [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ] = PolySolve_complex4( w(q),Ka,Ma,La,nor,1e-3);

    %3-1
    tol = 1e-5;
    lim = 2;
    lim2 = 2000;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]  = PolySolve_complex4( w(q),Ka,Ma,La,nor,tol,lim,lim2);

%     [ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex( w,Ka,Ma,La,nor,tol);
end

ai = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN

[a] = sortDiff(ai,f);
% a = ai;

for q=1:length(f)
    %% B
    %30-10
%     tol = 1e-4;
%     [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]  = PolySolve_complex4( w(q),Kb,Mb,Lb,nor,1e-3); %for 30-10 cases, use lim1 and lim2 default

    %3-1
    tol = 1e-5;
    lim = 0.1;
    lim2 = 10;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]  = PolySolve_complex4( w(q),Kb,Mb,Lb,nor,tol,lim,lim2);

%     [ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex( w,Kb,Mb,Lb,nor,tol);
end

bi = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
[b] = sortDiff(bi,f);
% b = bi;

for q=1:length(f)
    %% C
    %30-10
%     tol = 1e-4;
%     [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]  = PolySolve_complex4( w(q),Kc,Mc,Lc,nor,1e-3);
    %3-1
    tol = 1e-5;
    lim = 2;
    lim2 = 2000;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]  = PolySolve_complex4( w(q),Kc,Mc,Lc,nor,tol,lim,lim2);
    
%     [ PhiQc_p,PhiQc_n,PhiFc_p,PhiFc_n,PsiQc_p,PsiFc_p,PsiQc_n,PsiFc_n,lpc1,lnc1,sc,kpci,knci ] = PolySolve_complex( w,Kc,Mc,Lc,nor,tol);
    

end

ci = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
[c] = sortDiff(ci,f);    
% c = ci;

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
    
%     mode_limitb = 3;
%     mode_limitb = length(knbi);
    mode_limitb = nmodes_b;
    b.PhiQp = b.PhiQp(:,1:mode_limitb,:);
    b.PhiQn = b.PhiQn(:,1:mode_limitb,:);
    
    b.PsiQp = b.PsiQp(1:mode_limitb,:,:);
    b.PsiQn = b.PsiQn(1:mode_limitb,:,:);
    
    b.PhiFp = b.PhiFp(:,1:mode_limitb,:);
    b.PhiFn = b.PhiFn(:,1:mode_limitb,:);
    
%     PsiFb_p1(:,:,q) = PsiFb_p(1:mode_limitb,:);
%     PsiFb_n1(:,:,q) = PsiFb_n(1:mode_limitb,:);

    mode_limitc = nmodes_a;
    c.PhiQp = c.PhiQp(:,1:mode_limitc,:);
    c.PhiQn = c.PhiQn(:,1:mode_limitc,:);
    
    c.PsiQp = c.PsiQp(1:mode_limitc,:,:);
    c.PsiQn = c.PsiQn(1:mode_limitc,:,:);
    
    c.PhiFp = c.PhiFp(:,1:mode_limitc,:);
    c.PhiFn = c.PhiFn(:,1:mode_limitc,:);


for q=1:length(f)
    %% Analytical Solution
    
    kaa(q) = sqrt(w(q))*(rho*Sa/E/Ia)^(1/4);   % Wave number
    kbb(q) = sqrt(w(q))*(rho*Sb/E/Ib)^(1/4);   % Wave number
    kbc(q) = sqrt(w(q))*(rho*Sc/E/Ic)^(1/4); 
    
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
    
    %% Boundary conditions when using the left eigenvectors
    [ar,ac] = size(Ma);
    [br,bc] = size(Mb);
    
    Eb = zeros(ar/2,br/2);
    Ca = zeros(ar/2);
    
    % Describe the boundary (connection a-b)
    I = eye(bc/2);
    Inoda = (ar/2-br/2)/2+1:(ar/2-br/2)/2+br/2; %SM
%     Inoda = (ar/2-44/2)/2+1:(ar/2-44/2)/2+br/2; %Quasi SM
%     Inoda = 1:br/2; 
    Inodb = 1:br/2;
    
    Ea = eye(ar/2);
    Eb(Inoda,Inodb) = I;
    Ec = Ea;
    Ca(Inoda,Inoda) = I;
    Cb = Eb;
    Cc = Ca;
    
    md = 1;
    
    CP1ba = [PsiQa_p0*Ca*a.PhiQn(:,:,q) -PsiQa_p0*Cb*b.PhiQp(:,:,q);
        PsiFa_p0*Ea*a.PhiFn(:,:,q) -PsiFa_p0*Eb*b.PhiFp(:,:,q)];
    
    CP2ba = [-PsiQa_p0*Ca*a.PhiQp(:,:,q) PsiQa_p0*Cb*b.PhiQn(:,:,q);
        -PsiFa_p0*Ea*a.PhiFp(:,:,q) PsiFa_p0*Eb*b.PhiFn(:,:,q)];
    
%     CP1ba = [PsiQa_p0*Ca*a.PhiQn(:,:,q) -PsiQa_p0*Cb*b.PhiQp(:,:,q);
%         PsiFa_p0*Ea*a.PhiFn(:,:,q)  PsiFa_p0*Eb*b.PhiFp(:,:,q)];
%     
%     CP2ba = [-PsiQa_p0*Ca*a.PhiQp(:,:,q) PsiQa_p0*Cb*b.PhiQn(:,:,q);
%         -PsiFa_p0*Ea*a.PhiFp(:,:,q) -PsiFa_p0*Eb*b.PhiFn(:,:,q)];
    
    TRTba =(pinv(CP1ba)*CP2ba);
%     TRTba =(pinv(CP2ba)*CP1ba);
    [TRTr(q), TRTc(q)] = size(TRTba);
%     TRT2Pba(:,:,q) = TRTba;
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

%     CP1cb = [PsiQa_p0*Cb*PhiQb_n1 -PsiQa_p0*Cc*PhiQc_p1;
%         PsiFa_p0*Eb*PhiFb_n1 -PsiFa_p0*Ec*PhiFc_p1];
%     
%     CP2cb = [-PsiQa_p0*Cb*PhiQb_p1 PsiQa_p0*Cc*PhiQc_n1;
%         -PsiFa_p0*Eb*PhiFb_p1 PsiFa_p0*Ec*PhiFc_n1];
%     
    CP1cb = [PsiQa_p0*Cb*b.PhiQn(:,:,q) -PsiQa_p0*Cc*c.PhiQp(:,:,q);
        PsiFa_p0*Eb*b.PhiFn(:,:,q) -PsiFa_p0*Ec*c.PhiFp(:,:,q)];
    
    CP2cb = [-PsiQa_p0*Cb*b.PhiQp(:,:,q) PsiQa_p0*Cc*c.PhiQn(:,:,q);
        -PsiFa_p0*Eb*b.PhiFp(:,:,q) PsiFa_p0*Ec*c.PhiFn(:,:,q)];

%     CP1cb = [PsiQa_p0*Cb*b.PhiQn(:,:,q) -PsiQa_p0*Cc*c.PhiQp(:,:,q);
%         PsiFa_p0*Eb*b.PhiFn(:,:,q)  PsiFa_p0*Ec*c.PhiFp(:,:,q)];
%     
%     CP2cb = [-PsiQa_p0*Cb*b.PhiQp(:,:,q) PsiQa_p0*Cc*c.PhiQn(:,:,q);
%         -PsiFa_p0*Eb*b.PhiFp(:,:,q) -PsiFa_p0*Ec*c.PhiFn(:,:,q)];

%     TRTPcb(:,:,q) = (pinv(CP1cb)*CP2cb);
    TRTPcb = (pinv(CP1cb)*CP2cb);
%     TRTPcb = (pinv(CP2cb)*CP1cb);
    
    RPbb2(:,:,q) = TRTPcb(1:nmodes_b,1:nmodes_b);
    TPcb2(:,:,q) = TRTPcb(nmodes_b+1:end,1:nmodes_b);
    RPcc2(:,:,q) = TRTPcb(nmodes_b+1:end,nmodes_b+1:end);
    TPbc2(:,:,q) = TRTPcb(1:nmodes_b,nmodes_b+1:end);
    
    % Scattering from a to c
    [RPAA(:,:,q),TPCA(:,:,q)] = ThreeSectionRT(RPaa1(:,:,q),RPbb2(:,:,q),RPbb1(:,:,q),TPba1(:,:,q),TPcb2(:,:,q),TPab1(:,:,q),TPb);

    
    %%
    
     mode_limita = nmodes_a;
    PhiQa_p2 = a.PhiQp(:,1:mode_limita,q);
    PhiQa_n2 = a.PhiQn(:,1:mode_limita,q);
    
    PhiFa_p2 = a.PhiFp(:,1:mode_limita,q);
    PhiFa_n2 = a.PhiFn(:,1:mode_limita,q);
    
    lpa2 = a.lp(1:mode_limita,q);
    lna2 = a.ln(1:mode_limita,q);
    
%     mode_limitb = 3;
%     mode_limitb = length(knbi);
    mode_limitb = nmodes_b;
    PhiQb_p2 = b.PhiQp(:,1:mode_limitb,q);
    PhiQb_n2 = b.PhiQn(:,1:mode_limitb,q);
    
    PsiQb_p2(:,:,q) = b.PsiQp(1:mode_limitb,:,q);
    PsiQb_n2(:,:,q) = b.PsiQn(1:mode_limitb,:,q);
    
    PhiFb_p2 = b.PhiFp(:,1:mode_limitb,q);
    PhiFb_n2 = b.PhiFn(:,1:mode_limitb,q);
        
    lpb2 = b.lp(1:mode_limitb,q);
    lnb2 = b.ln(1:mode_limitb,q);
    
    %     mode_limitc = 3;
%     mode_limitc = length(knci);
    mode_limitc = nmodes_c;
    PhiQc_p2 = c.PhiQp(:,1:mode_limitc,q);
    PhiQc_n2 = c.PhiQn(:,1:mode_limitc,q);
    
    PsiQc_p2(:,:,q) = c.PsiQp(1:mode_limitc,:,q);
    PsiQc_n2(:,:,q) = c.PsiQn(1:mode_limitc,:,q);
    
    PhiFc_p2 = c.PhiFp(:,1:mode_limitc,q);
    PhiFc_n2 = c.PhiFn(:,1:mode_limitc,q);
        
    lpc2 = c.lp(1:mode_limitc,q);
    lnc2 = c.ln(1:mode_limitc,q);
    
    
     %% Power matrix for the Plane element 1
    % MITROU (2015)
    Pa2 = (1i*w(q)/2)*[PhiQa_p2'*PhiFa_p2 PhiQa_p2'*PhiFa_n2;
        PhiQa_n2'*PhiFa_p2 PhiQa_n2'*PhiFa_n2]-...
        [PhiFa_p2'*PhiQa_p2 PhiFa_p2'*PhiQa_n2;
        PhiFa_n2'*PhiQa_p2 PhiFa_n2'*PhiQa_n2];
    Pc2 = (1i*w(q)/2)*[PhiQc_p2'*PhiFc_p2 PhiQc_p2'*PhiFc_n2;
        PhiQc_n2'*PhiFc_p2 PhiQc_n2'*PhiFc_n2]-...
        [PhiFc_p2'*PhiQc_p2 PhiFc_p2'*PhiQc_n2;
        PhiFc_n2'*PhiQc_p2 PhiFc_n2'*PhiQc_n2];
    
%     PrPP2(q) = abs(TRT(1,1))^2*(Pa2(28,28)/Pa2(1,1));
%     PrPL2(q) = abs(TRT(2,1))^2*(Pa2(29,29)/Pa2(1,1));
%     PrPN2(q) = abs(TRT(3,1))^2*(Pa2(30,30)/Pa2(1,1));
%     PrPC12(q) = abs(TRT(4,1))^2*(Pa2(31,31)/Pa2(1,1));
%     PrPC22(q) = abs(TRT(5,1))^2*(Pa2(32,32)/Pa2(1,1));
%     PtPP2(q) = abs(TRT(nmodes_a+1,1))^2*(Pb2(1,1)/Pa2(1,1));
%     PtPL2(q) = abs(TRT(nmodes_a+2,1))^2*(Pb2(2,2)/Pa2(1,1));
%     PtPN2(q) = abs(TRT(nmodes_a+3,1))^2*(Pb2(3,3)/Pa2(1,1));
%     PtPC12(q) = abs(TRT(nmodes_a+4,1))^2*(Pb2(4,4)/Pa2(1,1));
%     PtPC22(q) = abs(TRT(nmodes_a+5,1))^2*(Pb2(5,5)/Pa2(1,1));
%    
%     PrLP2(q) = abs(TRT(1,2))^2*(Pa2(28,28)/Pa2(2,2));
%     PrLL2(q) = abs(TRT(2,2))^2*(Pa2(29,29)/Pa2(2,2));
%     PrLN2(q) = abs(TRT(3,2))^2*(Pa2(30,30)/Pa2(2,2));
%     PrLC12(q) = abs(TRT(4,2))^2*(Pa2(31,31)/Pa2(2,2));
%     PrLC22(q) = abs(TRT(5,2))^2*(Pa2(32,32)/Pa2(2,2));
%     PtLP2(q) = abs(TRT(nmodes_a+1,2))^2*(Pb2(1,1)/Pa2(2,2));
%     PtLL2(q) = abs(TRT(nmodes_a+2,2))^2*(Pb2(2,2)/Pa2(2,2));
%     PtLN2(q) = abs(TRT(nmodes_a+3,2))^2*(Pb2(3,3)/Pa2(2,2));
%     PtLC12(q) = abs(TRT(nmodes_a+4,2))^2*(Pb2(4,4)/Pa2(2,2));
%     PtLC22(q) = abs(TRT(nmodes_a+5,2))^2*(Pb2(5,5)/Pa2(2,2));
% 

    PrPP2(q) = abs(RPAA(1,1,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(1,1));
    PrPL2(q) = abs(RPAA(2,1,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(1,1));
    PrPN2(q) = abs(RPAA(3,1,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(1,1));
%     PrPC12(q) = abs(RPAA(4,1,q))^2*(Pa2(nmodes_a+4,nmodes_a+4)/Pa2(1,1));
%     PrPC22(q) = abs(RPAA(5,1,q))^2*(Pa2(nmodes_a+5,nmodes_a+5)/Pa2(1,1));
    PtPP2(q) = abs(TPCA(1,1,q))^2*(Pc2(1,1)/Pa2(1,1));
    PtPL2(q) = abs(TPCA(2,1,q))^2*(Pc2(2,2)/Pa2(1,1));
    PtPN2(q) = abs(TPCA(3,1,q))^2*(Pc2(3,3)/Pa2(1,1));
%     PtPC12(q) = abs(TPCA(4,1,q))^2*(Pc2(4,4)/Pa2(1,1));
%     PtPC22(q) = abs(TPCA(5,1,q))^2*(Pc2(5,5)/Pa2(1,1));
   
    PrLP2(q) = abs(RPAA(1,2,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(2,2));
    PrLL2(q) = abs(RPAA(2,2,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(2,2));
    PrLN2(q) = abs(RPAA(3,2,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(2,2));
%     PrLC12(q) = abs(RPAA(4,2,q))^2*(Pa2(nmodes_a+4,nmodes_a+4)/Pa2(2,2));
%     PrLC22(q) = abs(RPAA(5,2,q))^2*(Pa2(nmodes_a+5,nmodes_a+5)/Pa2(2,2));
    PtLP2(q) = abs(TPCA(1,2,q))^2*(Pc2(1,1)/Pa2(2,2));
    PtLL2(q) = abs(TPCA(2,2,q))^2*(Pc2(2,2)/Pa2(2,2));
    PtLN2(q) = abs(TPCA(3,2,q))^2*(Pc2(3,3)/Pa2(2,2));
%     PtLC12(q) = abs(TPCA(4,2,q))^2*(Pc2(4,4)/Pa2(2,2));
%     PtLC22(q) = abs(TPCA(5,2,q))^2*(Pc2(5,5)/Pa2(2,2));
    
    %% Kinetic Energy - BENDING INCIDENT
    
    AMP_I = zeros(nmodes_a,1);
    AMP_I(1) = 1;
    AMP_R = RPAA(:,:,q)*AMP_I;
%     AMP_T = zeros(nmodes_b,1);
%     AMP_T(1) = 1;
    AMP_T = TPCA(:,:,q)*AMP_I;
    
    % Bending Propagating Positive in A  (incident)
    [EIbbkin_avL(q),EIbbpot_av(q),EIbbpot_perx(q),EIbbpot_pery(q),PIbbp(q),PRbbx(q),PIbby(q),EIbbkin_avLx(q),EIbbkin_avLy(q),veIbbp(q)] = KineticEn3(a.PhiQp(:,:,q),a.lp(q),Ma,Ka,w(q),La,AMP_I);
    
    % Bending Propagating Negative in A  (reflected)
    [ERbbkin_avL(q),ERbbpot_av(q),ERbbpot_perx(q),ERbbpot_pery(q),PRbbn(q),PRbbx(q),PRbby(q),ERbbkin_avLx(q),ERbbkin_avLy(q),veRbbn(q)] = KineticEn3(a.PhiQn(:,:,q),a.ln(q),Ma,Ka,w(q),La,AMP_R);

    % Bending Propagating Positive in B (transmitted) 
    [ETbbkin_avL(q),ETbbpot_av(q),ETbbpot_perx(q),ETbbpot_pery(q),PTbbp(q),PTbbx(q),PTbby(q),ETbbkin_avLx(q),ETbbkin_avLy(q),veTbbp(q)]  = KineticEn3(c.PhiQp(:,:,q),c.lp(q),Mc,Kc,w(q),Lc,AMP_T);

    %% Kinetic Energy - LONGITUDINAL INCIDENT 
 
    AMP_I = zeros(nmodes_a,1);
    AMP_I(2) = 1;
    AMP_R = RPAA(:,:,q)*AMP_I;
%     AMP_T = zeros(nmodes_b,1);
%     AMP_T(2) = 1;
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
% R_WFE2(3,:) = reshape(TRT22(1,1,:),[1 length(f)]);
T_WFE2(1,:) = reshape(TPCA(1,1,:),[1 length(f)]);
T_WFE2(2,:) = reshape(TPCA(2,2,:),[1 length(f)]);
% T_WFE2(3,:) = reshape(TRT22(2,1,:),[1 length(f)]);

R_WFE1(1,:) = reshape(RPaa1(1,1,:),[1 length(f)]);
R_WFE1(2,:) = reshape(RPaa1(2,2,:),[1 length(f)]);
% R_WFE2(3,:) = reshape(TRT22(1,1,:),[1 length(f)]);
T_WFE1(1,:) = reshape(TPab1(1,1,:),[1 length(f)]);
T_WFE1(2,:) = reshape(TPba1(2,2,:),[1 length(f)]);

% Cutoff frequencies
Cutoff_ai = [88e3 154e3 166e3 178e3 266e3]; 
Cutoff_bi = 170e3;

Cutoff_a = round(Cutoff_ai/df+1); 
Cutoff_b = round(Cutoff_bi/df);
% Cutoff_a = round(([156e3 166e3 ]-fi)/df); 
% Cutoff_b = round((266e3-fi)/df);


%% Standard Plots
% 
% % This plot shows the reflection with longitudinal incident and
% % longitudinal reflected. B: analytical, R: numerical
%  figure()
% %  Analytical bending R
%  plot(f,abs(R_WM(1,:)),'b-','LineWidth',1)
%  hold on
% %  PLANE ELEMENT bending R
%  plot(f,abs(R_WFE2(1,:)),'b*','LineWidth',1)
% 
% %  Analytical long R
%  plot(f,abs(R_WM(2,:)),'r-.','LineWidth',1)
% 
% %  PLANE ELEMENT long R
%  plot(f,abs(R_WFE2(2,:)),'ro','LineWidth',1)
% 
% plot([5e3 5e3],[-800 800],'g:','LineWidth',1.5)
% plot([30e3 30e3],[-800 800],'g:','LineWidth',1.5)
% plot([54.5e3 54.5e3],[-800 800],'g:','LineWidth',1.5)
% plot([70e3 70e3],[-800 800],'g:','LineWidth',1.5)
% plot([88.5e3 88.5e3],[-800 800],'g:','LineWidth',1.5)
% plot([110e3 110e3],[-800 800],'g:','LineWidth',1.5)
% plot([150e3 150e3],[-800 800],'g:','LineWidth',1.5)
% plot([160e3 160e3],[-800 800],'g:','LineWidth',1.5)
% plot([166.5e3 166.5e3],[-800 800],'g:','LineWidth',1.5)
% plot([171.5e3 171.5e3],[-800 800],'g:','LineWidth',1.5)
% plot([177.5e3 177.5e3],[-800 800],'g:','LineWidth',1.5)
% plot([250e3 250e3],[-800 800],'g:','LineWidth',1.5)
% plot([300e3 300e3],[-800 800],'g:','LineWidth',1.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff 0 3])
% 
%  figure()
% % Analytical bending T
% plot(f,abs(T_WM(1,:)),'b-','LineWidth',1)
% hold on
% % Plane bending T
% plot(f,abs(T_WFE2(1,:)),'b*','LineWidth',1)
% 
% % Analytical long T
% plot(f,abs(T_WM(2,:)),'r-.','LineWidth',1)
% 
% % Plane long T
% plot(f,abs(T_WFE2(2,:)),'ro','LineWidth',1)
% 
% plot([5e3 5e3],[-3 3],'g:','LineWidth',1.5)
% plot([30e3 30e3],[-3 3],'g:','LineWidth',1.5)
% plot([54.5e3 54.5e3],[-3 3],'g:','LineWidth',1.5)
% plot([70e3 70e3],[-3 3],'g:','LineWidth',1.5)
% plot([88.5e3 88.5e3],[-3 3],'g:','LineWidth',1.5)
% plot([110e3 110e3],[-3 3],'g:','LineWidth',1.5)
% plot([150e3 150e3],[-3 3],'g:','LineWidth',1.5)
% plot([160e3 160e3],[-3 3],'g:','LineWidth',1.5)
% plot([166.5e3 166.5e3],[-3 3],'g:','LineWidth',1.5)
% plot([171.5e3 171.5e3],[-3 3],'g:','LineWidth',1.5)
% plot([177.5e3 177.5e3],[-3 3],'g:','LineWidth',1.5)
% plot([250e3 250e3],[-3 3],'g:','LineWidth',1.5)
% plot([300e3 300e3],[-3 3],'g:','LineWidth',1.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff 0 4])
% 
% 
% 
%  figure()
% %  Analytical bending R
%  plot(f,phase(R_WM(1,:))*180/pi,'b-','LineWidth',1)
%  hold on
% %  PLANE ELEMENT bending R
%  plot(f,atan2(imag(R_WFE2(1,:)),real(R_WFE2(1,:)))*180/pi,'b*','LineWidth',1)
% 
%  %  Analytical long R
%  plot(f,phase(R_WM(2,:))*180/pi,'r-.','LineWidth',1)
% 
% %  PLANE ELEMENT long R
%  plot(f,atan2(imag(R_WFE2(2,:)),real(R_WFE2(2,:)))*180/pi,'ro','LineWidth',1)
% 
% % legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','\phi_R [deg]','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff -180 180])
% %
%  figure()
%  % Analytical bending T
%  plot(f,phase(T_WM(1,:))*180/pi,'b-','LineWidth',1)
%  hold on
%  % Plane bending T
%  plot(f,atan2(imag(T_WFE2(1,:)),real(T_WFE2(1,:)))*180/pi,'b*','LineWidth',1)
% 
%  % Analytical long T
%  plot(f,phase(T_WM(2,:))*180/pi,'r-.','LineWidth',1)
% 
%  % Plane long T
%  plot(f,atan2(imag(T_WFE2(2,:)),real(T_WFE2(2,:)))*180/pi,'ro','LineWidth',1)
% 
% % legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','\phi_T [deg]','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff -180 180])

%% NON Cutoff plots
% % Plot this first to check which coefficients are zero. This helps to have
% % less lines in the plot
% figure(1)
% plot(f,abs(reshape(RPAA(1,1,:),[1 length(f)])));
% hold on
% plot(f,abs(reshape(RPAA(2,1,:),[1 length(f)])));
% plot(f,abs(reshape(RPAA(3,1,:),[1 length(f)])));
% plot(f,abs(reshape(RPAA(4,1,:),[1 length(f)])));
% plot(f,abs(reshape(RPAA(5,1,:),[1 length(f)])));
% 
% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([70e3 70e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff 0 3])
%  
% figure(2)
% plot(f,abs(reshape(TPCA(1,1,:),[1 length(f)])));
% hold on
% plot(f,abs(reshape(TPCA(2,1,:),[1 length(f)])));
% plot(f,abs(reshape(TPCA(3,1,:),[1 length(f)])));
% plot(f,abs(reshape(TPCA(4,1,:),[1 length(f)])));
% plot(f,abs(reshape(TPCA(5,1,:),[1 length(f)])));
% 
% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([70e3 70e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff 0 3])
%  
%  figure(3)
% plot(f,abs(reshape(RPAA(1,2,:),[1 length(f)])));
% hold on
% plot(f,abs(reshape(RPAA(2,2,:),[1 length(f)])));
% plot(f,abs(reshape(RPAA(3,2,:),[1 length(f)])));
% plot(f,abs(reshape(RPAA(4,2,:),[1 length(f)])));
% plot(f,abs(reshape(RPAA(5,2,:),[1 length(f)])));
% 
% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([70e3 70e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff 0 3])
%  
%  figure(4)
% plot(f,abs(reshape(TPCA(1,2,:),[1 length(f)])));
% hold on
% plot(f,abs(reshape(TPCA(2,2,:),[1 length(f)])));
% plot(f,abs(reshape(TPCA(3,2,:),[1 length(f)])));
% plot(f,abs(reshape(TPCA(4,2,:),[1 length(f)])));
% plot(f,abs(reshape(TPCA(5,2,:),[1 length(f)])));
% 
% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([70e3 70e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([fi ff 0 3])
%% Cutoff plots
% Use the plotCutOffs only for the non zero coefficients, the others plot
% them normally

figure()

plot(f,abs(reshape(RPAA(1,1,:),[1 length(f)])),'b','LineWidth',2);
hold on
% plot(f,abs(reshape(RPAA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(RPAA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(RPAA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPAA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);
% plotCutOffs(f,abs(reshape(RPAA(6,1,:),[1 length(f)])),Cutoff_a(5), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPAA(7,1,:),[1 length(f)])),Cutoff_a(5), {'g','r-.'},[1,2]);

plot(f,abs(reshape(RBTAA(1,1,:),[1 length(f)])),'b:')
% plot(f,abs(reshape(RBTAA(2,1,:),[1 length(f)])),'b:')

for ii = 1:length(Cutoff_ai)
    plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
end
for ii = 1:length(Cutoff_bi)
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

% plot(f,abs(reshape(TPCA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(TPCA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(TPCA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPCA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);
% plotCutOffs(f,abs(reshape(TPCA(6,1,:),[1 length(f)])),Cutoff_a(5), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPCA(7,1,:),[1 length(f)])),Cutoff_a(5), {'g','r-.'},[1,2]);

plot(f,abs(reshape(TBTCA(1,1,:),[1 length(f)])),'b:')
% plot(f,abs(reshape(TBTCA(2,1,:),[1 length(f)])),'b:')

for ii = 1:length(Cutoff_ai)
    plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
end
for ii = 1:length(Cutoff_bi)
    plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
end

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
% plotCutOffs(f,abs(reshape(RPAA(3,2,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(RPAA(4,2,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPAA(5,2,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);
% plotCutOffs(f,abs(reshape(RPAA(6,1,:),[1 length(f)])),Cutoff_a(5), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPAA(7,1,:),[1 length(f)])),Cutoff_a(5), {'g','r-.'},[1,2]);
% 
% plot(f,abs(R_WM(2,:,:)),'b:')
% 
% for ii = 1:length(Cutoff_ai)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_bi)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
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
% plotCutOffs(f,abs(reshape(TPCA(6,1,:),[1 length(f)])),Cutoff_a(5), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPCA(7,1,:),[1 length(f)])),Cutoff_a(5), {'g','r-.'},[1,2]);
% 
% plot(f,abs(T_WM(2,:)),'b:')
% 
% for ii = 1:length(Cutoff_ai)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_bi)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])
 
 %% IND
%  
% figure()
% hold on
% plot(f,abs(reshape(RPaa1(1,1,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(RPaa1(2,1,:),[1 length(f)])),'b','LineWidth',2);
% % plotCutOffs(f,abs(reshape(RPaa1(3,1,:),[1 length(f)])),Cutoff_a(1), {'r-.','b'},[2,2]);
% % plotCutOffs(f,abs(reshape(RPbb2(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% % plotCutOffs(f,abs(reshape(RPbb2(5,1,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);
% 
% plot(f,abs(reshape(RBTaa1(1,1,:),[1 length(f)])),'b:')
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
%  axis([fi ff 0 3])
%  
%  figure()
% hold on
% plot(f,abs(reshape(TPba1(1,1,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(TPba1(2,1,:),[1 length(f)])),'b','LineWidth',2);
% % plotCutOffs(f,abs(reshape(TPba1(3,1,:),[1 length(f)])),Cutoff_a(1), {'r-.','b'},[2,2]);
% % plotCutOffs(f,abs(reshape(TPbc2(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% % plotCutOffs(f,abs(reshape(TPbc2(5,1,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);
% 
% plot(f,abs(reshape(TBTba1(1,1,:),[1 length(f)])),'b:')
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
%  axis([fi ff 0 3])
% 
%  figure()
% hold on
% plot(f,abs(reshape(RPaa1(1,2,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(RPaa1(2,2,:),[1 length(f)])),'b','LineWidth',2);
% % plotCutOffs(f,abs(reshape(RPaa1(3,2,:),[1 length(f)])),Cutoff_a(1), {'r-.','b'},[2,2]);
% % plotCutOffs(f,abs(reshape(RPbb2(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% % plotCutOffs(f,abs(reshape(RPbb2(5,1,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);
% 
% plot(f,abs(reshape(RLTaa1,[1 length(f)])),'b:')
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
%  axis([fi ff 0 3])
%  
%  figure()
% hold on
% plot(f,abs(reshape(TPba1(1,2,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(TPba1(2,2,:),[1 length(f)])),'b','LineWidth',2);
% % plotCutOffs(f,abs(reshape(TPba1(3,2,:),[1 length(f)])),Cutoff_a(1), {'r-.','b'},[2,2]);
% % plotCutOffs(f,abs(reshape(TPbc2(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% % plotCutOffs(f,abs(reshape(TPbc2(5,1,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);
% 
% plot(f,abs(reshape(TLTba1,[1 length(f)])),'b:')
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
%  axis([fi ff 0 3])
 %% 
%  figure()
% hold on
% plot(f,abs(reshape(RPaa1(1,1,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(RPaa1(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(RPaa1(3,1,:),[1 length(f)])),Cutoff_a(1), {'r-.','b'},[2,2]);
% plot(f,abs(reshape(RPaa1(4,1,:),[1 length(f)])),'g','LineWidth',1);
% plot(f,abs(reshape(RPaa1(5,1,:),[1 length(f)])),'g','LineWidth',1);
% 
% plot(f,abs(reshape(RBTaa1(1,1,:),[1 length(f)])),'b:')
% 
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
%  axis([fi ff 0 3])
% 
% figure()
% hold on
% plot(f,abs(reshape(TPba1(1,1,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(TPba1(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(TPba1(3,1,:),[1 length(f)])),Cutoff_a(1), {'r-.','b'},[2,2]);
% plot(f,abs(reshape(TPba1(4,1,:),[1 length(f)])),'g','LineWidth',1);
% plot(f,abs(reshape(TPba1(5,1,:),[1 length(f)])),'g','LineWidth',1);
% 
% plot(f,abs(reshape(TBTba1(1,1,:),[1 length(f)])),'b:')
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
%  axis([fi ff 0 3])
%  
%  figure()
% plot(f,abs(reshape(RPaa1(1,2,:),[1 length(f)])),'b','LineWidth',2);
% hold on
% plot(f,abs(reshape(RPaa1(2,2,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(RPaa1(3,2,:),[1 length(f)])),'r-.','LineWidth',2);
% plotCutOffs(f,abs(reshape(RPaa1(4,2,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPaa1(5,2,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);
% 
% plot(f,abs(RLTaa1(1,:)),'b:')
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
%  axis([fi ff 0 3])
% 
% figure()
% plot(f,abs(reshape(TPba1(1,2,:),[1 length(f)])),'b','LineWidth',2);
% hold on
% plot(f,abs(reshape(TPba1(2,2,:),[1 length(f)])),'b','LineWidth',2);
% plot(f,abs(reshape(TPba1(3,2,:),[1 length(f)])),'r-.','LineWidth',2);
% plotCutOffs(f,abs(reshape(TPba1(4,2,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPba1(5,2,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);
% 
% n = length(f);
% plot(f,abs(TLTba1(1,:)),'b:')
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
%  axis([fi ff 0 3])
%  
%% Energy plots
% 
%  figure()
%  plot(f,abs(EIbbkin_avLy./EIbbkin_avL)*100,'b*','LineWidth',1)
%  hold on
%  plot(f,abs(EIbbkin_avLx./EIbbkin_avL)*100,'rp','LineWidth',1)
%  legend('kEy_{inc}/kE_{inc}','kEx_{inc}/kE_{inc}')
% %  title('Symmetric discontinuity')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Time-averaged kinetic energy [%]','FontName','Times New Roman','FontSize',12)
%  set(gca,'fontsize',12);
%  axis([fi f(q) -1 101])
%  
%  for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
% % Bending Incident Power Coefficients
%  figure()
%  plot(f,abs(PrPP2),'b','LineWidth',1)
%  hold on
%  plot(f,abs(PrPL2),'r','LineWidth',1)
%  plot(f,abs(PrPN2),'g','LineWidth',1)
% %  plot(f,abs(PrPC12),'m','LineWidth',1)
% %  plot(f,abs(PrPC22),'c','LineWidth',1)
%  plot(f,abs(PtPP2),'b*','LineWidth',1)
%  plot(f,abs(PtPL2),'r*','LineWidth',1)
%  plot(f,abs(PtPN2),'g*','LineWidth',1)
% %  plot(f,abs(PtPC12),'m*','LineWidth',1)
% %  plot(f,abs(PtPC22),'c*','LineWidth',1)
% 
% for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
%  legend('Propagating bending reflected','Longitudinal reflected','Nearfield reflected','Attenuated 1 reflected','Attenuated 2 reflected','Propagating bending transmitted','Longitudinal transmitted','Nearfield transmitted','Attenuated 1 transmitted','Attenuated 2 transmitted')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Power coefficients','FontName','Times New Roman','FontSize',12)
%  set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi f(q) -0.01 2])
% 
%  % Bending Incident TA Kinetic Energy
%  figure()
%  plot(f,abs(EIbbkin_avLy./EIbbkin_avL)*100,'b*','LineWidth',1)
%  hold on
%  plot(f,abs(EIbbkin_avLx./EIbbkin_avL)*100,'rp','LineWidth',1)
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
%  plot(f,abs(ERbbkin_avLy./ERbbkin_avL)*100,'b*','LineWidth',1)
%  hold on
%  plot(f,abs(ERbbkin_avLx./ERbbkin_avL)*100,'rp','LineWidth',1)
%  plot(f,abs(ETbbkin_avLy./ETbbkin_avL)*100,'go','LineWidth',1)
%  plot(f,abs(ETbbkin_avLx./ETbbkin_avL)*100,'ms','LineWidth',1)
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
% % Longitudinal Incident Power Coefficients
%  figure()
%  plot(f,abs(PrLP2),'b','LineWidth',1)
%  hold on
%  plot(f,abs(PrLL2),'r','LineWidth',1)
%  plot(f,abs(PrLN2),'g','LineWidth',1)
% %  plot(f,abs(PrLC12),'m','LineWidth',1)
% %  plot(f,abs(PrLC22),'c','LineWidth',1)
%  plot(f,abs(PtLP2),'b*','LineWidth',1)
%  plot(f,abs(PtLL2),'r*','LineWidth',1)
%  plot(f,abs(PtLN2),'g*','LineWidth',1)
% %  plot(f,abs(PtLC12),'m*','LineWidth',1)
% %  plot(f,abs(PtLC22),'c*','LineWidth',1)
%  for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
%  legend('Propagating bending reflected','Longitudinal reflected','Nearfield reflected','Attenuated 1 reflected','Attenuated 2 reflected','Propagating bending transmitted','Longitudinal transmitted','Nearfield transmitted','Attenuated 1 transmitted','Attenuated 2 transmitted')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Power coefficients','FontName','Times New Roman','FontSize',12)
%  set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi f(q) -0.01 2])
% 
%   % Longitudinal Incident TA Kinetic Energy
%  figure()
%  plot(f,abs(EIllkin_avLy./EIllkin_avL)*100,'b*','LineWidth',1)
%  hold on
%  plot(f,abs(EIllkin_avLx./EIllkin_avL)*100,'rp','LineWidth',1)
% 
%  for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
%  legend('kEy_{inc}/kE_{inc}','kEx_{inc}/kE_{inc}')
% % %  title('Symmetric discontinuity')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Time-averaged kinetic energy [%]','FontName','Times New Roman','FontSize',12)
%  set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi f(q) -1 101])
% 
%  figure()
%  plot(f,abs(ERlbkin_avLy./ERlbkin_avL)*100,'b*','LineWidth',1)
%  hold on
%  plot(f,abs(ERlbkin_avLx./ERlbkin_avL)*100,'rp','LineWidth',1)
%  plot(f,abs(ETlbkin_avLy./ETlbkin_avL)*100,'go','LineWidth',1)
%  plot(f,abs(ETlbkin_avLx./ETlbkin_avL)*100,'ms','LineWidth',1)
%  
%  for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% 
%  legend('kEy_{refl}/kE_{refl}','kEx_{refl}/kE_{refl}','kEy_{trans}/kE_{trans}','kEx_{trans}/kE_{trans}')
% % %  title('Symmetric discontinuity')
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Time-averaged kinetic energy [%]','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi f(q) -1 101])
%  
%% dispersion plot
% figure()
% plot(f,a.kp(:,:),'b')
% hold on
% % plot (f, kbaf(:),'r')
% % plot (f, klaf(:),'r')
% plot(f,imag(a.kp(:,:)),'b')
% xlabel('Frequency (Hz)');
% ylabel('Real(k)');
% set(gca,'fontsize',12,'FontName','Times New Roman');
% 
% % plot(1:10)
% 
% box off
% 
% axes('xlim', [0 3e5], 'ylim', [-700 700], 'color', 'none', 'YAxisLocation', 'right', 'XAxisLocation', 'top','YTick',[],'XTick',[])
% set(gca,'fontsize',12,'FontName','Times New Roman');
% ylabel('Imag(k)');
% 
% % 
% 
% figure()
% plot(f,b.kp(:,:),'b')
% hold on
% % plot (f, kbaf(:),'r')
% % plot (f, klaf(:),'r')
% plot(f,imag(b.kp(:,:)),'b')
% xlabel('Frequency (Hz)');
% ylabel('Real(k)');
% set(gca,'fontsize',12,'FontName','Times New Roman');
% 
% % plot(1:10)
% 
% box off
% 
% axes('xlim', [0 3e5], 'ylim', [-700 700], 'color', 'none', 'YAxisLocation', 'right', 'XAxisLocation', 'top','YTick',[],'XTick',[])
% set(gca,'fontsize',12,'FontName','Times New Roman');
% ylabel('Imag(k)');
% % 
% % 


%% Simple dispersion
% figure()
% plot(f,a.kp(:,:),'b')
% xlabel('Frequency (Hz)');
% ylabel('Real(k)');
% set(gca,'fontsize',12,'FontName','Times New Roman');
% 
% figure()
% plot(f,imag(a.kp(:,:)),'b')
% xlabel('Frequency (Hz)');
% ylabel('Imag(k)');
% set(gca,'fontsize',12,'FontName','Times New Roman');
% 
% figure()
% plot(f,b.kp(:,:),'b')
% xlabel('Frequency (Hz)');
% ylabel('Real(k)');
% set(gca,'fontsize',12,'FontName','Times New Roman');
% 
% figure()
% plot(f,imag(b.kp(:,:)),'b')
% xlabel('Frequency (Hz)');
% ylabel('Imag(k)');
% set(gca,'fontsize',12,'FontName','Times New Roman');

%% Colored dispersion

figure()
plot(f,real(a.kp(1:2,:)),'b','LineWidth',2);
hold on
plotCutOffs(f,real(a.kp(3,:)),Cutoff_a(1), {'r-.','b'},[2,2]);
plotCutOffs(f,imag(a.kp(3,:)),Cutoff_a(1), {'r-.','b'},[2,2]);
plotCutOffs(f,real(a.kp(4,:)),Cutoff_a(2), {'g','b'},[1,2]);
plotCutOffs(f,imag(a.kp(4,:)),Cutoff_a(2), {'m','b'},[1,2]);
plotCutOffs(f,real(a.kp(5,:)),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);
plotCutOffs(f,imag(a.kp(5,:)),Cutoff_a(2:4), {'m','b','r-.','b'},[1,2,2,2]);
plotCutOffs(f,real(a.kp(6,:)),Cutoff_a(5), {'g','b'},[1,2]);
plotCutOffs(f,imag(a.kp(6,:)),Cutoff_a(5), {'m','r-.'},[1,2]);
plotCutOffs(f,real(a.kp(7,:)),Cutoff_a(5), {'g','b'},[1,2]);
plotCutOffs(f,imag(a.kp(7,:)),Cutoff_a(5), {'m','r-.'},[1,2]);

plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
plot([0 ff],[0 0],'k','LineWidth',1)

xlabel('Frequency (Hz)');
ylabel('Real(k)');
set(gca,'fontsize',12,'FontName','Times New Roman');

figure()
plot(f,real(b.kp(1:2,:)),'b','LineWidth',2);
hold on
plotCutOffs(f,real(b.kp(3,:)),Cutoff_a(5), {'r-.','b'},[2,2]);
plotCutOffs(f,imag(b.kp(3,:)),Cutoff_a(5), {'r-.','b'},[2,2]);
plot(f,real(b.kp(4:5,:)),'g','LineWidth',1);
plot(f,imag(b.kp(4:5,:)),'m','LineWidth',1);

plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)
plot([0 ff],[0 0],'k','LineWidth',1)

xlabel('Frequency (Hz)');
ylabel('Real(k)');
set(gca,'fontsize',12,'FontName','Times New Roman');
