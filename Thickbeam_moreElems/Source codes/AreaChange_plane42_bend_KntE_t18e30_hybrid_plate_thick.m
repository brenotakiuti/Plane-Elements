% Find the reflection and transmission from two sections modeled in FE
% using commercial package
% long
% Breno Ebinuma Takiuti
% 03/11/2017

clear 
clc
% close all


%% Beam Analytical solution

%% Material's Constants
% Material: Steel

rho=7800;     %mass per unit valume
E=2.06e11;  %Young's modulus

%% Geometric constants

Ndof = 2;
b = 6e-3;              % base of the cross-section (tickness) (m)
ha = 18e-3;              % height of the cross-section (width) (m)
hb = 6e-3;
% hb = 18e-3;
% hb = 36e-3;
hc = 18e-3;
NeleB = 150;            % Number of elements in B
La = 6e-4;
Lb = 6e-4;               % length of the element (x direction) (m) (3 elem)
Lc = La;
L = NeleB*Lb;         % Use NeleB to calculate a length of B
Sa = b*(ha);
Sb = b*hb;
Sc = b*hc;
Ia = b*(ha)^3/12;
Ib = b*(hb)^3/12;
Ic = b*(hc)^3/12;
beta_ab = Sb/Sa;
beta_bc = Sc/Sb;

nModes = 2;

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

Kabc = zeros(ra+(ra/2*(NeleB+1)));
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
% ra
for ii=1:NeleB
    Inodab(:,:,ii) = [round(ra/2*ii)+Inoda, ra+(ra/2*(ii-1))+Inoda];
end
Inodbc = [ra+1+(ra/2*(NeleB-1)):ra+(ra/2*(NeleB+1))];

Ea = eye(ar/2);
Eb(Inoda,Inodb) = I;
Ec = Ea;
Ca(Inoda,Inoda) = I;
Cb = Eb;
Cc = Ca;

md = 1;

for ii=1:NeleB
    Kabc( Inodab(:,:,ii), Inodab(:,:,ii)) = Kabc( Inodab(:,:,ii), Inodab(:,:,ii))+Kb;
    Mabc( Inodab(:,:,ii), Inodab(:,:,ii)) = Mabc( Inodab(:,:,ii), Inodab(:,:,ii))+Mb;
end
Kabc(Inodbc,Inodbc) = Kabc(Inodbc,Inodbc)+Kc0;
Mabc(Inodbc,Inodbc) = Mabc(Inodbc,Inodbc)+Mc0;

%%
% Frequencies
fi = 100;
ff = 250000;
df = 1e4;
f = fi:df:ff;
% f = 10e3;

% w = f(1)*2*pi;

%% Normalization
nor = 0;

%% FE matrices

% [Ma,Ka] = EB_Beam(rho,Sa,La,E,Ia);
% [Mb,Kb] = EB_Beam(rho,Sb,Lb,E,Ib);
% [Mb2,Kb2] = EB_Beam(rho,Sb,L,E,Ib);
% [Mc,Kc] = EB_Beam(rho,Sc,Lc,E,Ic);
    
%%

tol = 1e-4;

w = 2*pi*f;

for q=1:length(f)

    %% Waveguide A
    tol = 1e-5;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
      = PolySolve_complex6( w(q),Ka,Ma,La,nor,tol);

% [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
%       PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
%       = PolySolve_complex6( w(q),Ka,Ma,La,nor,tol);
end

ai = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN

tolS = 1e-1;
% [a] = sortDiff(ai,f,tolS);
a = ai;
%%
for q=1:length(f)
    %% B
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
    tol = 1e-5;
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
    nmodes_a = length(a.kp(:,1));
    nmodes_b = 3;
    nmodes_c = length(c.kp(:,1));
    
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
    
    AA = eye(2);
    BB = AA;

    kaa(q) = sqrt(w(q))*(rho*Sa/E/Ia)^(1/4);   % Wave number
    kbb(q) = sqrt(w(q))*(rho*Sb/E/Ib)^(1/4);   % Wave number
    kbc(q) = sqrt(w(q))*(rho*Sc/E/Ic)^(1/4); 
    
%     lbpaa=exp(-1i*kaa*La);
%     lbnaa=exp(1i*kaa*La);
%     lbpab=exp(-1i*kbb*Lb);
%     lbpac=exp(-1i*kbc*Lc);
%     
    
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
end

for q=1:length(f)
    %% Numeric
            
%     nmodes_b = length(b.kp(:,q));
%     nmodes_a = length(a.kp(:,q));
%     kPb = b.kp (1:nmodes_b,q);
    
%     PsiQa_p0 = a.PsiQp(:,:,1);
%     PsiFa_p0 = a.PsiFp(:,:,1);
    

    %% Hybrid method by Renno
    
    Z1 = zeros(size(a.PhiQp(:,:,q)));
    Z2 = Z1';
    I = eye(size(a.PhiQp(:,:,q)));
    [rp,cp,~] = size(a.PhiQp);
    [rb,cb,~] = size(Kb);
  
    PhiQ_p = [a.PhiQp(:,:,q) Z1; Z1 c.PhiQp(:,:,q)];
    PhiQ_n = [a.PhiQn(:,:,q) Z1; Z1 c.PhiQn(:,:,q)];
    PhiF_p = [a.PhiFp(:,:,q) Z1; Z1 c.PhiFp(:,:,q)];
    PhiF_n = [a.PhiFn(:,:,q) Z1; Z1 c.PhiFn(:,:,q)];
    

    PsiQ_n = [a.PsiQn(:,:,q) Z2; Z2 c.PsiQn(:,:,q)];
%     PsiQ_n = [I Z2; Z2 I];

%     RA = -eye(rp); RB = -eye(rp);
    RA = eye(rp); RB = eye(rp);
%     RA = eye(rp); RB = -eye(rp);
    Z3 = zeros(rp);
    
    R = [RA Z3; Z3 RB];
    
    Nind = rp+1:2:rp*2;
    
    R(Nind,Nind) = -eye(length(Nind));
    

    Kabc2 = zeros(NeleB*Ndof+2);
    Mabc2 = Kabc2;

    D0 = (Kabc2-w(q)^2*Mabc2)*1;
    D1 = (Kabc-w(q)^2*Mabc)*1;
    D2 = (Ka-w(q)^2*Ma);
    
%    [ra,ca] = size(Kabc2);
   [rabc, cabc] = size(Kabc);
    
    InodE = [1:ra/2,rabc-ra/2+1:rabc];
    InodI = ra/2+1:rabc-ra/2;

%     %Method 1
    DEE = D1(InodE,InodE);
    DEI = D1(InodE,InodI);
    DIE = D1(InodI,InodE);
    DII = D1(InodI,InodI);
    
    Dj = (DEE-DEI*pinv(DII)*DIE);    %Method 1

    %Method 2
%     DTLL = D0(InodE(1:2),InodE(1:2));
%     DTLR = D0(InodE(1:2),InodE(3:4));
%     DTRL = D0(InodE(3:4),InodE(1:2));
%     DTRR = D0(InodE(3:4),InodE(3:4));
%     DTLO = D0(InodE(1:2),InodI);
%     DTRO = D0(InodE(3:4),InodI);
%     DTOL = D0(InodI,InodE(1:2));
%     DTOR = D0(InodI,InodE(3:4));
%     DTOO = D0(InodI,InodI);
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

    
end
[rp,cp,~] = size(a.PhiQp);
RHAA(:,:,:) = Srt(1:cp,1:cp,:);
THCA(:,:,:) = Srt(cp+1:cp+nModes,1:cp,:);
% TRTP = [RPAA; TPCA];

% Analyticlal Reflection and Transmission
% R_WM(1,:) = RBTAA(1,1,:);
% T_WM(1,:) = TBTCA(1,1,:);
% R_WM(2,:) = RLTAA;
% T_WM(2,:) = TLTCA;
% R_WM(3,:) = RBTAA(2,1,:);
% T_WM(3,:) = TBTCA(2,1,:);

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
% plot(f,abs(reshape(RPAA(1,1,:),[1 length(f)])),'b','LineWidth',2);
% hold on
% plot(f,abs(reshape(RPAA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(RPAA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(RPAA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(RPAA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);

plot(f,abs(reshape(RBTAA(1,1,:),[1 length(f)])),'b:')
% plot(f,abs(reshape(RBTAA(2,1,:),[1 length(f)])),'b:')
hold on
plot(f,abs(reshape(RHAA(1,1,:),[1 length(f)])),'r--')
% plot(f,abs(reshape(RHAA(1,2,:),[1 length(f)])),'r--')
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
% plot(f,abs(reshape(TPCA(1,1,:),[1 length(f)])),'b','LineWidth',2);
% hold on
% plot(f,abs(reshape(TPCA(2,1,:),[1 length(f)])),'b','LineWidth',2);
% plotCutOffs(f,abs(reshape(TPCA(3,1,:),[1 length(f)])),88, {'r-.','b'},[2,2]);
% plotCutOffs(f,abs(reshape(TPCA(4,1,:),[1 length(f)])),Cutoff_a(2), {'g','b'},[1,2]);
% plotCutOffs(f,abs(reshape(TPCA(5,1,:),[1 length(f)])),Cutoff_a(2:4), {'g','b','r-.','b'},[1,2,2,2]);

plot(f,abs(reshape(TBTCA(1,1,:),[1 length(f)])),'b:')
% plot(f,abs(reshape(TBTCA(2,1,:),[1 length(f)])),'b:')
hold on
plot(f,abs(reshape(THCA(1,1,:),[1 length(f)])),'r--')
% plot(f,abs(reshape(THCA(1,2,:),[1 length(f)])),'r--')
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
