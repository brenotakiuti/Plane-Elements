% Find the reflection and transmission from two sections modeled in FE
% using commercial package
% bend: both long and bending
% non dimensional plot
% Breno Ebinuma Takiuti
% 23/08/2017

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
Lb = 6e-4;               % length of the element (m)
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
% load KM_plane42_Lx001Ly0075_12el.mat  %12 Plane elements with Lx = 0.001 Ly=0.0075

% load KM_plane43_Lx0006Ly006_10el.mat
% load KM_plane43_Lx0006Ly012_20el.mat
load KM_plane43_Lx0006Ly018_30el.mat
Ka=full(K);
Ma=full(M);

Kc = Ka;
Mc = Ma;

% K and M for the second section
% load KM_plane42_Lx001Ly0025_4el.mat  %4 Plane elements with Lx=0.001 Ly=0.0025

load KM_plane43_Lx0006Ly006_10el.mat
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
Kb=full(K);
Mb=full(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2; %degree of freedom per node
na=30; %number of structural elements for the first section
nb=20; %number of elements for the second section
nc=na;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dofa = n*(na+1);       % number of degrees of freedom per side
dofb = n*(nb+1);
dofc = dofa;

%% Normalization
nor = 0;

%%

% tol = 1e-4; %S3 and NS3 for 30-10
tol = 1e-3;
w0 = 100*2*pi;    

[ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex( w0,Ka,Ma,La,nor,tol);
nmodes_a = length(kpai);
% mode_limita = 5;
PsiQa_p0 = PsiQa_p(1:nmodes_a,:);
PsiQa_n0 = PsiQa_n(1:nmodes_a,:);
PsiFa_p0 = PsiFa_p(1:nmodes_a,:);
PsiFa_n0 = PsiFa_n(1:nmodes_a,:);

[ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex( w0,Kb,Mb,Lb,nor,tol);
% nmodes_b = length(kpbi)
nmodes_b = 3;       % Observation later
PsiQb_p0 = PsiQb_p(1:nmodes_b,:);
PsiQb_n0 = PsiQb_n(1:nmodes_b,:);
PsiFb_p0 = PsiFb_p(1:nmodes_b,:);
PsiFb_n0 = PsiFb_n(1:nmodes_b,:);

nmodes_c = nmodes_a;

% Frequencies
fi = 100;
ff = 300000;
df = 100;
xx = sqrt(fi):sqrt(ff);
f = xx.^2;
% f = fi:df:ff;
endf = 20e3;


w = 2*pi*f;
wf = 2*pi*endf;

ka = sqrt(w)*(rho*Sa/E/Ia)^(1/4);   % Wave number
kaf = sqrt(wf)*(rho*Sa/E/Ia)^(1/4);   % Wave number

% kaf = sqrt( 2*pi*endf)*(rho*Sa/E/Ia)^(1/4);

lambda_l = 2*pi./ka;
lambda_lf = 2*pi./kaf;

% lambda_lf = 2*pi./sqrt( 2*pi*f)*(rho*Sa/E/Ia)^(1/4);
L = min(lambda_lf);
% L = 0.018;
% L = 1e-6;

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

    %% Waveguide A

%     [ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex( w(q),Ka,Ma,La,nor,tol);
    [ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,lpa1,lna1,sa,kpai,knai ] = PolySolve_complex( w(q),Ka,Ma,La,nor,1e-3,4,650);
    kpa(:,q) = kpai;
    kna(:,q) = knai;
    sar(:,q) = sa;
    
    ppa(q) = length(kpai);
    
    %% B
%     [ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex( w(q),Kb,Mb,Lb,nor,tol); 
    [ PhiQb_p,PhiQb_n,PhiFb_p,PhiFb_n,PsiQb_p,PsiFb_p,PsiQb_n,PsiFb_n,lpb1,lnb1,sb,kpbi,knbi ] = PolySolve_complex( w(q),Kb,Mb,Lb,nor,1e-3,4,600); %for 30-10 cases, use lim1 and lim2 default
    kpb(:,q) = kpbi;
    knb(:,q) = knbi;
    kPb = 1i*log(lpb1)/Lb;
    
    ppb(q) = length(kpbi);
    
    %% C
%     [ PhiQc_p,PhiQc_n,PhiFc_p,PhiFc_n,PsiQc_p,PsiFc_p,PsiQc_n,PsiFc_n,lpc1,lnc1,sc,kpci,knci ] = PolySolve_complex( w(q),Kc,Mc,Lc,nor,tol);
    [ PhiQc_p,PhiQc_n,PhiFc_p,PhiFc_n,PsiQc_p,PsiFc_p,PsiQc_n,PsiFc_n,lpc1,lnc1,sc,kpci,knci ] = PolySolve_complex( w(q),Kc,Mc,Lc,nor,1e-3,4,650);
  
    kpc(:,q) = kpci;
    knc(:,q) = knci;
    
    ppc(q) = length(kpci);
    
    %% Should always use all propagating and nearfield waves, cant cut
%     mode_limita = 3;
%     
%     PhiQa_p1 = PhiQa_p(:,1:nmodes_a);
%     PhiQa_n1 = PhiQa_n(:,1:nmodes_a);
%     
%     PhiFa_p1 = PhiFa_p(:,1:nmodes_a);
%     PhiFa_n1 = PhiFa_n(:,1:nmodes_a);
%     
%     mode_limitb = 3;
%     mode_limitb = length(knbi);
%     PhiQb_p1 = PhiQb_p(:,1:mode_limitb);
%     PhiQb_n1 = PhiQb_n(:,1:mode_limitb);
%     
%     PsiQb_p1(:,:,q) = PsiQb_p(1:mode_limitb,:);
%     PsiQb_n1(:,:,q) = PsiQb_n(1:mode_limitb,:);
%     
%     PhiFb_p1 = PhiFb_p(:,1:mode_limitb);
%     PhiFb_n1 = PhiFb_n(:,1:mode_limitb);
%     
% %     PsiFb_p1(:,:,q) = PsiFb_p(1:mode_limitb,:);
% %     PsiFb_n1(:,:,q) = PsiFb_n(1:mode_limitb,:);
%     

%     mode_limita = 3;
    mode_limita = nmodes_a;
    PhiQa_p1 = PhiQa_p(:,1:mode_limita);
    PhiQa_n1 = PhiQa_n(:,1:mode_limita);
    
    PhiFa_p1 = PhiFa_p(:,1:mode_limita);
    PhiFa_n1 = PhiFa_n(:,1:mode_limita);
    
    lpa1 = lpa1(1:mode_limita);
    lna1 = lna1(1:mode_limita);
    
%     mode_limitb = 3;
%     mode_limitb = length(knbi);
    mode_limitb = nmodes_b;
    PhiQb_p1 = PhiQb_p(:,1:mode_limitb);
    PhiQb_n1 = PhiQb_n(:,1:mode_limitb);
    
    PsiQb_p1(:,:,q) = PsiQb_p(1:mode_limitb,:);
    PsiQb_n1(:,:,q) = PsiQb_n(1:mode_limitb,:);
    
    PhiFb_p1 = PhiFb_p(:,1:mode_limitb);
    PhiFb_n1 = PhiFb_n(:,1:mode_limitb);
        
    lpb1 = lpb1(1:mode_limitb);
    lnb1 = lnb1(1:mode_limitb);
    
%         mode_limitc = 3;
%     mode_limitc = length(knci);
    mode_limitc = nmodes_c;
    PhiQc_p1 = PhiQc_p(:,1:mode_limitc);
    PhiQc_n1 = PhiQc_n(:,1:mode_limitc);
    
    PsiQc_p1(:,:,q) = PsiQc_p(1:mode_limitc,:);
    PsiQc_n1(:,:,q) = PsiQc_n(1:mode_limitc,:);
    
    PhiFc_p1 = PhiFc_p(:,1:mode_limitc);
    PhiFc_n1 = PhiFc_n(:,1:mode_limitc);
        
    lpc1 = lpc1(1:mode_limitc);
    lnc1 = lnc1(1:mode_limitc);

%     nmodes_a = mode_limita;
%     nmodes_b = mode_limitb;
%     nmodes_c = mode_limitc;
%     
    kPb = kPb (1:nmodes_b);
    %% Boundary conditions when using the left eigenvectors
    [ar,ac] = size(Ma);
    [br,bc] = size(Mb);
    
    Eb = zeros(ar/2,br/2);
    Ca = zeros(ar/2);
    
    % Describe the boundary (connection a-b)
    I = eye(bc/2);
%     Inoda = (ar/2-br/2)/2+1:(ar/2-br/2)/2+br/2; %SM
    Inoda = 1:br/2; 
    Inodb = 1:br/2;
    
    Ea = eye(ar/2);
    Eb(Inoda,Inodb) = I;
    Ec = Ea;
    Ca(Inoda,Inoda) = I;
    Cb = Eb;
    Cc = Ca;
    
    md = 1;
    
%     C12 = [PsiQa_p*Ca*PhiQa_n1 -PsiQa_p*Cb*PhiQb_p1;
%         PsiFa_p*Ea*PhiFa_n1 -PsiFa_p*Eb*PhiFb_p1];
%     
%     C22 = [-PsiQa_p*Ca*PhiQa_p1 PsiQa_p*Cb*PhiQb_n1;
%         -PsiFa_p*Ea*PhiFa_p1 PsiFa_p*Eb*PhiFb_n1];
    
    CP1ba = [PsiQa_p0*Ca*PhiQa_n1 -PsiQa_p0*Cb*PhiQb_p1;
        PsiFa_p0*Ea*PhiFa_n1 -PsiFa_p0*Eb*PhiFb_p1];
    
    CP2ba = [-PsiQa_p0*Ca*PhiQa_p1 PsiQa_p0*Cb*PhiQb_n1;
        -PsiFa_p0*Ea*PhiFa_p1 PsiFa_p0*Eb*PhiFb_n1];
    %
    %     C12 = [Ca*PhiQa_n1 -Cb*PhiQb_p1;
    %         Ea*PhiFa_n1 -Eb*PhiFb_p1];
    %
    %     C22 = [-Ca*PhiQa_p1 Cb*PhiQb_n1;
    %         -Ea*PhiFa_p1 Eb*PhiFb_n1];
    
    TRTba =(pinv(CP1ba)*CP2ba);
    [TRTr(q), TRTc(q)] = size(TRTba);
    TRT2Pba(:,:,q) = TRTba;
    
    mode = round(TRTc(q)/2);
    
    RPaa1(:,:,q) = TRT2Pba(1:nmodes_a,1:nmodes_a,q);
    TPba1(:,:,q) = TRT2Pba(nmodes_a+1:end,1:nmodes_a,q);
    RPbb1(:,:,q) = TRT2Pba(nmodes_a+1:end,nmodes_a+1:end,q);
    TPab1(:,:,q) = TRT2Pba(1:nmodes_a,nmodes_a+1:end,q);
    
    % Bending transition matrix from interface 1 to 2
    TPb = diag(exp(-1i*kPb*L));
    
    %% Important observation:
    % At first we where using 5 modes in waveguide B, 2 of them being
    % complex modes. I found that using these modes to construct the TPb
    % matrix generates numerical errors. If we cut to 3 modes (propagating
    % bending, longitudinal and nearfield), the results does not change 
    % much as compared to cutting the number of modes in A, while avoiding 
    % this error.
    
    %%
    
    % Matrix inversion method as in Harland et al (2000) eq. (54) - b to c
%     CP1cb = [PsiQa_p0*Cb*PhiQb_n1 -PsiQa_p0*Cc*PhiQc_p1;
%         PsiFa_p0*Eb*PhiFb_n1 -PsiFa_p0*Ec*PhiFc_p1];
%     
%     CP2cb = [-PsiQa_p0*Cb*PhiQb_p1 PsiQa_p0*Cc*PhiQc_n1;
%         -PsiFa_p0*Eb*PhiFb_p1 PsiFa_p0*Ec*PhiFc_n1];

    CP1cb = [PsiQa_p0*Cb*PhiQb_n1 -PsiQa_p0*Cc*PhiQc_p1;
        PsiFa_p0*Eb*PhiFb_n1 -PsiFa_p0*Ec*PhiFc_p1];
    
    CP2cb = [-PsiQa_p0*Cb*PhiQb_p1 PsiQa_p0*Cc*PhiQc_n1;
        -PsiFa_p0*Eb*PhiFb_p1 PsiFa_p0*Ec*PhiFc_n1];
    
    TRTPcb(:,:,q) = (pinv(CP1cb)*CP2cb);
    
    RPbb2(:,:,q) = TRTPcb(1:nmodes_b,1:nmodes_b,q);
    TPcb2(:,:,q) = TRTPcb(nmodes_b+1:end,1:nmodes_b,q);
    RPcc2(:,:,q) = TRTPcb(nmodes_b+1:end,nmodes_b+1:end,q);
    TPbc2(:,:,q) = TRTPcb(1:nmodes_b,nmodes_b+1:end,q);
    
    % Scattering from a to c
    [RPAA(:,:,q),TPCA(:,:,q)] = ThreeSectionRT(RPaa1(:,:,q),RPbb2(:,:,q),RPbb1(:,:,q),TPba1(:,:,q),TPcb2(:,:,q),TPab1(:,:,q),TPb);
%     [dofb,~]=size(RPbb1(:,:,q));
%     
%     RPAA(:,:,q) = RPaa1(:,:,q) + TPab1(:,:,q)*pinv(eye(dofb)-TPb*RPbb2(:,:,q)*TPb*RPbb1(:,:,q))*TPb*RPbb2(:,:,q)*TPb*TPba1(:,:,q);
%     TPCA(:,:,q) = TPcb2(:,:,q)*pinv(eye(dofb)-TPb*RPbb1(:,:,q)*TPb*RPbb2(:,:,q))*TPb*TPba1(:,:,q);
  
    TRTP = [RPAA; TPCA];

end

% mode_limita = 5;

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

% Cutoff frequencies
% Cutoff_a = [round(88000/df) round(156000/df) round(166000/df) round(178000/df)]; 
% Cutoff_b = round(266000/df);
Cutoff_a = [88e3 156e3 166e3 178e3];
Cutoff_b = 266e3;

Cutoff_a = round(sqrt(Cutoff_a))-xx(1)+1;
Cutoff_b = round(sqrt(Cutoff_b))-xx(1)+1;

%% Standard Plots

% This plot shows the reflection with longitudinal incident and
% longitudinal reflected. B: analytical, R: numerical

%  figure()
% %  Analytical bending R
%  plot(L./lambda_l,abs(R_WM(1,:)),'b-','LineWidth',1)
%  hold on
% %  PLANE ELEMENT bending R
%  plot(L./lambda_l,abs(R_WFE2(1,:)),'r-.','LineWidth',1)
%  % Analytical bending T
% plot(L./lambda_l,abs(T_WM(1,:)),'b-','LineWidth',1)
% % Plane bending T
% plot(L./lambda_l,abs(T_WFE2(1,:)),'r-.','LineWidth',1)
% 
% % Black lines
% plot([0 max(L./lambda_l)],[1 1],'k:')
% 
% legend('Analytical results', 'Plane elements results')
%  set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Reflection and transmission coefficients','FontName','Times New Roman','FontSize',12)
% set(gca,'FontName','Times New Roman','fontsize',12);
%  axis([0 max(L./lambda_l) 0 1.2])
%  box off
%  ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  axis([fi ff 0 1.2])
%  set(gca,'fontsize',12,'FontName','Times New Roman');
%  
%   figure()
% %  Analytical bending R
%  plot(L./lambda_l,atan2(imag(R_WM(1,:)),real(R_WM(1,:)))*180/pi,'b-','LineWidth',1)
%  hold on
% %  PLANE ELEMENT bending R
%  plot(L./lambda_l,atan2(imag(R_WFE2(1,:)),real(R_WFE2(1,:)))*180/pi,'r-.','LineWidth',1)
% %  % Analytical bending T
%  plot(L./lambda_l,atan2(imag(T_WM(1,:)),real(T_WM(1,:)))*180/pi,'b-','LineWidth',1)
% %  % Plane bending T
%  plot(L./lambda_l,atan2(imag(T_WFE2(1,:)),real(T_WFE2(1,:)))*180/pi,'r-.','LineWidth',1)
% % plot(L./lambda_l,phase(T_WFE2(1,:))*180/pi,'r-.','LineWidth',1)
% 
% % Black lines
% plot([0 max(L./lambda_l)],[0 0],'k-','LineWidth',0.8)
% 
% legend('Analytical results', 'Plane elements results')
%  set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','Phase [deg]','FontName','Times New Roman','FontSize',12)
% set(gca,'FontName','Times New Roman','fontsize',12);
%  axis([0 max(L./lambda_l) -270 270])
%  box off
%  ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  axis([fi ff  -270 270])
%  set(gca,'fontsize',12,'FontName','Times New Roman');
 
%  figure()
% %  Analytical bending R
%  plot(L./lambda_l,abs(R_WM(1,:)),'b-','LineWidth',1)
%  hold on
% %  PLANE ELEMENT bending R
%  plot(L./lambda_l,abs(R_WFE2(1,:)),'b*','LineWidth',1)
% 
% %  Analytical long R
%  plot(L./lambda_l,abs(R_WM(2,:)),'r-.','LineWidth',1)
% 
% %  PLANE ELEMENT long R
%  plot(L./lambda_l,abs(R_WFE2(2,:)),'ro','LineWidth',1)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([0 max(L./lambda_l) 0 2])
%  box off
%  ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
%   set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
% 
%  axis([fi ff 0 1])
%  set(gca,'fontsize',12,'FontName','Times New Roman');
%  
%  figure()
% % Analytical bending T
% plot(L./lambda_l,abs(T_WM(1,:)),'b-','LineWidth',1)
% hold on
% % Plane bending T
% plot(L./lambda_l,abs(T_WFE2(1,:)),'b*','LineWidth',1)
% 
% % Analytical long T
% plot(L./lambda_l,abs(T_WM(2,:)),'r-.','LineWidth',1)
% 
% % Plane long T
% plot(L./lambda_l,abs(T_WFE2(2,:)),'ro','LineWidth',1)
% 
% %  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
% %  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([0 max(L./lambda_l) 0 2])
%  box off
%  ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
%  set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
%  axis([fi ff 0 1])
%  set(gca,'fontsize',12,'FontName','Times New Roman');

% 
%  figure()
% %  Analytical bending R
%  plot(L./lambda_l,phase(R_WM(1,:))*180/pi,'b-','LineWidth',1)
%  hold on
% %  PLANE ELEMENT bending R
%  plot(L./lambda_l,unwrap(atan2(imag(R_WFE2(1,:)),real(R_WFE2(1,:)))*180/pi),'b*','LineWidth',1)
% 
%  %  Analytical long R
%  plot(L./lambda_l,phase(R_WM(2,:))*180/pi,'r-.','LineWidth',1)
% 
% %  PLANE ELEMENT long R
%  plot(L./lambda_l,unwrap(atan2(imag(R_WFE2(2,:)),real(R_WFE2(2,:)))*180/pi),'ro','LineWidth',1)
% 
% % legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','\phi_R [deg]','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([0 max(L./lambda_l) -180 180])
% % %
%  figure()
%  % Analytical bending T
%  plot(L./lambda_l,phase(T_WM(1,:))*180/pi,'b-','LineWidth',1)
%  hold on
%  % Plane bending T
%  plot(L./lambda_l,unwrap(atan2(imag(T_WFE2(1,:)),real(T_WFE2(1,:)))*180/pi),'b*','LineWidth',1)
% 
%  % Analytical long T
%  plot(L./lambda_l,phase(T_WM(2,:))*180/pi,'r-.','LineWidth',1)
% 
%  % Plane long T
%  plot(L./lambda_l,unwrap(atan2(imag(T_WFE2(2,:)),real(T_WFE2(2,:)))*180/pi),'ro','LineWidth',1)
% 
% % legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
%  set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
%  set(get(gca,'YLabel'),'String','\phi_T [deg]','FontName','Times New Roman','FontSize',12)
% set(gca,'fontsize',12);
%  axis([0 max(L./lambda_l) -180 180])

%% NON Cutoff plots
% Plot this first to check which coefficients are zero. This helps to have
% less lines in the plot
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
% 

figure()
hold on
plot(L./lambda_l,abs(reshape(RPAA(1,1,:),[1 length(f)])),'b','LineWidth',2);
plot(L./lambda_l,abs(reshape(RPAA(2,1,:),[1 length(f)])),'b','LineWidth',2);
plotCutOffs(L./lambda_l,abs(reshape(RPAA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
plotCutOffs(L./lambda_l,abs(reshape(RPAA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
plotCutOffs(L./lambda_l,abs(reshape(RPAA(5,1,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);

plot(L./lambda_l,abs(reshape(RBTAA(1,1,:),[1 length(f)])),'b:')
plot(L./lambda_l,abs(reshape(RBTAA(2,1,:),[1 length(f)])),'b:')

% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
 axis([0 max(L./lambda_l) 0 2])
 box off
 ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 axis([fi ff 0 2])
 set(gca,'fontsize',12,'FontName','Times New Roman');
 
figure()
hold on
plot(L./lambda_l,abs(reshape(TPCA(1,1,:),[1 length(f)])),'b','LineWidth',2);
plot(L./lambda_l,abs(reshape(TPCA(2,1,:),[1 length(f)])),'b','LineWidth',2);
plotCutOffs(L./lambda_l,abs(reshape(TPCA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
plotCutOffs(L./lambda_l,abs(reshape(TPCA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
plotCutOffs(L./lambda_l,abs(reshape(TPCA(5,1,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);

plot(L./lambda_l,abs(reshape(TBTCA(1,1,:),[1 length(f)])),'b:')
plot(L./lambda_l,abs(reshape(TBTCA(2,1,:),[1 length(f)])),'b:')

% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
 axis([0 max(L./lambda_l) 0 2])
 box off
 ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 axis([fi ff 0 2])
 set(gca,'fontsize',12,'FontName','Times New Roman');
 
 
 figure()
plot(L./lambda_l,abs(reshape(RPAA(1,2,:),[1 length(f)])),'b','LineWidth',2);
hold on
plot(L./lambda_l,abs(reshape(RPAA(2,2,:),[1 length(f)])),'b','LineWidth',2);
plotCutOffs(L./lambda_l,abs(reshape(RPAA(3,2,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
plotCutOffs(L./lambda_l,abs(reshape(RPAA(4,2,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
plotCutOffs(L./lambda_l,abs(reshape(RPAA(5,2,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);

plot(L./lambda_l,abs(reshape(RLTAA(1,1,:),[1 length(f)])),'b:')

% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
 axis([0 max(L./lambda_l) 0 2])
 box off
 ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 axis([fi ff 0 2])
 set(gca,'fontsize',12,'FontName','Times New Roman');
 
 
figure()
plot(L./lambda_l,abs(reshape(TPCA(1,2,:),[1 length(f)])),'b','LineWidth',2);
hold on
plot(L./lambda_l,abs(reshape(TPCA(2,2,:),[1 length(f)])),'b','LineWidth',2);
plotCutOffs(L./lambda_l,abs(reshape(TPCA(3,2,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
plotCutOffs(L./lambda_l,abs(reshape(TPCA(4,2,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
plotCutOffs(L./lambda_l,abs(reshape(TPCA(5,2,:),[1 length(f)])),Cutoff_a(2:end), {'g','b','r-.','b'},[1,2,2,2]);

plot(L./lambda_l,abs(reshape(TLTCA(1,1,:),[1 length(f)])),'b:')

% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','l_D/\lambda','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
 axis([0 max(L./lambda_l) 0 2])
  box off
 ax2 = axes('XAxisLocation','top', 'YAxisLocation', 'right','Color','none', 'YTick',[]);
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 axis([fi ff 0 2])
 set(gca,'fontsize',12,'FontName','Times New Roman');

%% dispersion plot
figure()
% plot(f,kpb(:,:),'b')
hold on
% plot (f, kbb(:),'r')
% plot (f, klb(:),'r')
plot(f,imag(kpb(:,:)),'b')
xlabel('Frequency (Hz)');
ylabel('Real(k)');
set(gca,'fontsize',12,'FontName','Times New Roman');

% plot(1:10)

box off

axes('xlim', [0 3e5], 'ylim', [-700 700], 'color', 'none', 'YAxisLocation', 'right', 'XAxisLocation', 'top','YTick',[],'XTick',[])
set(gca,'fontsize',12,'FontName','Times New Roman');
ylabel('Imag(k)');

figure()
plot(f,kpa(:,:),'b')
hold on
% plot (f, kbaf(:),'r')
% plot (f, klaf(:),'r')
plot(f,imag(kpa(:,:)),'b')
xlabel('Frequency (Hz)');
ylabel('Real(k)');
set(gca,'fontsize',12,'FontName','Times New Roman');

% plot(1:10)

box off

axes('xlim', [0 3e5], 'ylim', [-700 700], 'color', 'none', 'YAxisLocation', 'right', 'XAxisLocation', 'top','YTick',[],'XTick',[])
set(gca,'fontsize',12,'FontName','Times New Roman');
ylabel('Imag(k)'); 