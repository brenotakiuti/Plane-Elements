% Find the reflection and transmission from two sections modeled in FE
% using commercial package
% long
% Breno Ebinuma Takiuti
% 03/11/2017

clear 
% clc
% close all


%% Beam Analytical solution

%% Material's Constants
% Material: Steel

rho=7800;     %mass per unit valume
E=2.06e11;  %Young's modulus

%% Geometric constants

b = 6e-3;              % base of the cross-section (tickness) (m)
ha = 18e-3;              % height of the cross-section (width) (m)
% ha = 6e-3;
hb = 6e-3;
hc = 18e-3;
multi = 100;
L = multi*6e-4;
La = 6e-4;
Lb = L/multi;               % length of the element (x direction) (m) (3 elem)
Lc = La;
Sa = b*(ha);
Sb = b*hb;
Sc = b*hc;
Ia = b*(ha)^3/12;
Ib = b*(hb)^3/12;
Ic = b*(hc)^3/12;
beta_ab = Sb/Sa;
beta_bc = Sc/Sb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2; %degree of freedom per node
na=30; %number of structural elements for the first section
nb=10; %number of elements for the second section
nc=na;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dofa = n*(na+1);       % number of degrees of freedom per side
dofb = n*(nb+1);
dofc = dofa;

nModes = 1;

%%
% Frequencies
fi = 100;
ff = 100000;
df = 1e2;
f = fi:df:ff;
% f = 10e3;

% w = f(1)*2*pi;

%% Normalization
nor = 0;

%% FE matrices

Ka2 = E*Sa/La*[1 -1; -1 1];
Ma2 = rho*Sa*La*[2/6 1/6; 1/6 2/6];

Kb2 = E*Sb./Lb*[1 -1; -1 1];
Mb2 = rho*Sb.*Lb*[2/6 1/6; 1/6 2/6];

Kc2 = E*Sc./Lc*[1 -1; -1 1];
Mc2 = rho*Sc.*Lc*[2/6 1/6; 1/6 2/6];

Kabc = E*Sb/L*[1 -1; -1 1];
Mabc = rho*Sb.*L*[2/6 1/6; 1/6 2/6];
%%

tol = 1e-4;

w = 2*pi*f;

Lab = L+La;

for q=1:length(f)

    %% Waveguide A
    tol = 1e-2;
%     lim1 = 1;
%     lim2 = 1;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
      = PolySolve_complex6( w(q),Ka2,Ma2,La,nor,tol);
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
      = PolySolve_complex6( w(q),Kb2,Mb2,Lb,nor,tol); %for 30-10 cases, use lim1 and lim2 default
end

bi = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
% [b] = sortDiff(bi,f);
b = bi;

for q=1:length(f)
    %% C
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
      = PolySolve_complex6( w(q),Kc2,Mc2,Lc,nor,tol);

end

ci = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
% [c] = sortDiff(ci,f);    
c = ci;

for q=1:length(f)
    %% Analytical Solution

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
     
    nmodes_b = length(b.kp(:,q));
    nmodes_a = 1;
    kPb = b.kp (1:nmodes_b,q);
    
    PsiQa_p0 = 1;
    PsiFa_p0 = 1;
    
    Ca = 1;
    Cb = 1;
    Ea = 1;
    Eb = 1;
    
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

    Cc = 1;
    Ec = 1;
      
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
    
%     PsiQ_n = [a.PsiQn(:,:,1) Z2; Z2 c.PsiQn(:,:,1)];
    PsiQ_n = [1 Z2; Z2 1];

    RA = eye(rp); RB = eye(rp);
%     RA = eye(rp); RB = eye(rp);
    Z3 = zeros(rp);
    
    R = [RA Z3; Z3 RB];
    
%     Kabc2 = Kb2;
%     Mabc2 = Mb2;

    Kabc2 = zeros(multi+1);
    Mabc2 = Kabc2;
    
    for i=1:multi
        Kabc2(i:i+1,i:i+1) = Kabc2(i:i+1,i:i+1) + Kb2;       
        Mabc2(i:i+1,i:i+1) = Mabc2(i:i+1,i:i+1) + Mb2;
    end

% Kabc(1:2,1:2) = Kabc(1:2,1:2) + Kb2;
% Mabc(1:2,1:2) = Mabc(1:2,1:2) + Mb2;
% Kabc(2:3,2:3) = Kabc(2:3,2:3) + Kb2;
% Mabc(2:3,2:3) = Mabc(2:3,2:3) + Mb2;
% Kabc(3:4,3:4) = Kabc(3:4,3:4) + Kb2;
% Mabc(3:4,3:4) = Mabc(3:4,3:4) + Mb2;

%     D0 = (Kabc-w(q)^2*Mabc)*1000;
    D0 = (Kabc2-w(q)^2*Mabc2)*1;
    D1 = (Kabc-w(q)^2*Mabc)*1;
    D2 = (Ka2-w(q)^2*Ma2);
    
    %Method 1
    [ra,ca] = size(Kabc2);
    
    InodE = [1,ra];
    InodI = 2:ra-1;

    DEE = D0(InodE,InodE);
    DEI = D0(InodE,InodI);
    DIE = D0(InodI,InodE);
    DII = D0(InodI,InodI);
    
    Dj = (DEE-DEI*DII^-1*DIE);    %Method 1

    %Method 2
%     DTLL = D0(1:ra/2,1:ra/2);
%     DTLR = D0(1:ra/2,ra*3/2+1:ra*2);
%     DTRL = D0(ra*3/2+1:ra*2,1:rp);
%     DTRR = D0(ra*3/2+1:ra*2,ra*3/2+1:ra*2);
%     DTLO = D0(1:ra/2,ra/2+1:ra*3/2);
%     DTRO = D0(ra*3/2+1:ra*2,ra/2+1:ra*3/2);
%     DTOL = D0(ra/2+1:ra*3/2,1:ra/2);
%     DTOR = D0(ra/2+1:ra*3/2,ra*3/2+1:ra*2);
%     DTOO = D0(ra/2+1:ra*3/2,ra/2+1:ra*3/2);
%     
%     DLL = DTLL - DTLO*pinv(DTOO)*DTOL;
%     DLR = DTLR - DTLO*pinv(DTOO)*DTOR;
%     DRL = DTRL - DTRO*pinv(DTOO)*DTOL;
%     DRR = DTRR - DTRO*pinv(DTOO)*DTOR;
%     
%     Dj = [DLL DLR; DRL DRR];    %Method 2

    %Method 3
%     Dj = D0;
        
    Srt(:,:,q) = pinv(PsiQ_n*(-Dj*R*PhiQ_n+R*PhiF_n))*PsiQ_n*(Dj*R*PhiQ_p-R*PhiF_p);
%     Srt(:,:,q) = -pinv(Dj*R*PhiQ_n-R*PhiF_n)*(Dj*R*PhiQ_p-R*PhiF_p);

    RHAA(:,:,q) = Srt(1:cp,1:cp,q);
    THCA(:,:,q) = Srt(cp+1:cp+nModes,1:cp,q);
    

end

TRTP = [RPAA; TPCA];

% PLANE1 Extract the R&T from scaterring matrix
R_WFE2(1,:) = reshape(RPAA(1,1,:),[1 length(f)]);
% R_WFE2(2,:) = reshape(RPAA(2,2,:),[1 length(f)]);
T_WFE2(1,:) = reshape(TPCA(1,1,:),[1 length(f)]);
% T_WFE2(2,:) = reshape(TPCA(2,2,:),[1 length(f)]);

R_WFE1(1,:) = reshape(RPaa1(1,1,:),[1 length(f)]);
% R_WFE1(2,:) = reshape(RPaa1(2,2,:),[1 length(f)]);
T_WFE1(1,:) = reshape(TPab1(1,1,:),[1 length(f)]);
% T_WFE1(2,:) = reshape(TPba1(2,2,:),[1 length(f)]);

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
% plot(f,abs(reshape(RPAA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(RPAA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(RPAA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPAA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);

plot(f,abs(reshape(RLTAA(1,1,:),[1 length(f)])),'b:')
% plot(f,abs(reshape(RBTAA(2,1,:),[1 length(f)])),'b:')

plot(f,abs(reshape(RHAA(1,1,:),[1 length(f)])),'r--')
% plot(f,abs(reshape(RHAA(2,1,:),[1 length(f)])),'r--')

% for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

figure()
plot(f,abs(reshape(TPCA(1,1,:),[1 length(f)])),'b','LineWidth',2);
hold on
% plot(f,abs(reshape(TPCA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(TPCA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(TPCA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPCA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);

plot(f,abs(reshape(TLTCA(1,1,:),[1 length(f)])),'b:')
% plot(f,abs(reshape(TBTCA(2,1,:),[1 length(f)])),'b:')

plot(f,abs(reshape(THCA(1,1,:),[1 length(f)])),'r--')
% plot(f,abs(reshape(THCA(2,1,:),[1 length(f)])),'r--')

% plot([88e3 88e3],[-800 800],'k:','LineWidth',0.5)
% plot([156e3 156e3],[-800 800],'k:','LineWidth',0.5)
% plot([166e3 166e3],[-800 800],'k:','LineWidth',0.5)
% plot([178e3 178e3],[-800 800],'k:','LineWidth',0.5)
% plot([266e3 266e3],[-800 800],'k:','LineWidth',0.5)

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])
 