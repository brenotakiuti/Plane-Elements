% Modular Coding:
% Plane element , solve eig and save
% Breno Ebinuma Takiuti
% 14/12/2017

clear 
clc
% close all

%% Geometric constants

% La = 6e-4;
% Lb = 6e-4;               % length of the element (x direction) (m) (3 elem)
% Lc = La;
La = 5e-3;
Lb = 5e-3;               % length of the element (x direction) (m) (3 elem)
Lc = La;
% La = 1e-3;
% Lb = 1e-3;               % length of the element (x direction) (m) (3 elem)
% Lc = La;

%% mass and stiffness from Ansys - plane
% K and M for the first section
load KM_plane42_Lx005Ly0075_3el.mat
% load KM_plane42_Lx001Ly0075_12el.mat  %12 Plane elements with Lx = 0.001 Ly=0.0075
% load KM_plane43_Lx0006Ly009_15el.mat
% load KM_plane43_Lx0006Ly006_10el.mat
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
% load KM_plane43_Lx0006Ly006_10el.mat
Ka=full(K);
Ma=full(M);
Kc=Ka;
Mc=Ma;

% K and M for the second section
% load KM_plane42_Lx001Ly0025_4el.mat
load KM_plane42_Lx005Ly0025_1el.mat
% load KM_plane43_Lx0006Ly006_10el.mat            % Cutoff 266kHz
% load KM_plane43_Lx0006Ly0072_12el.mat         % Cutoff 222kHz
% load KM_plane43_Lx0006Ly009_15el.mat          % Cutoff 170kHz
% load KM_plane43_Lx0006Ly012_20el.mat
% load KM_plane43_Lx0006Ly018_30el.mat
Kb=full(K);
Mb=full(M);

[ndofa,ac] = size(Ma);
[ndofb,bc] = size(Mb);
[ndofc,cc] = size(Ma);

%% Frequencies
fi = 100;
ff = 1000;
df = 1e0;
f = fi:df:ff;
w = 2*pi*f;

%% Normalization
nor = 0;
tol = 1e-5;

for q=1:length(f)

    %% Waveguide A
    %30-10
%     tol = 1e-4;
%     [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
%       PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
%       = PolySolve_complex( w(q),Ka,Ma,La,nor,tol);

    %3-1
    tol = 1e-5;
    lim = 2;
    lim2 = 2000;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
      = PolySolve_complex( w(q),Ka,Ma,La,nor,tol,lim,lim2);

end
ndof = ndofa;
l = La;
ai = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN,ndof,l);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN

tolS = 1e-1;
% [a] = sortDiff(ai,f,tolS);
a = ai;
%%
for q=1:length(f)
    %% B
    %30-10
%     tol = 1e-4;
%     [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
%       PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ] ...
%       = PolySolve_complex( w(q),Kb,Mb,Lb,nor,tol); %for 30-10 cases, use lim1 and lim2 default

    %3-1
    tol = 1e-5;
    lim = 0.1;
    lim2 = 10;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ] ...
      = PolySolve_complex( w(q),Kb,Mb,Lb,nor,tol,lim,lim2); %for 30-10 cases, use lim1 and lim2 default

end
ndof = ndofb;
l = Lb;
bi = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN,ndof,l);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
% [b] = sortDiff(bi,f);
b = bi;

for q=1:length(f)
    %% C
    %30-10
%     tol = 1e-4;
%     [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
%       PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
%       = PolySolve_complex( w(q),Kc,Mc,Lc,nor,tol);   
  
    %3-1
    tol = 1e-5;
    lim = 2;
    lim2 = 2000;
    [ PhiQp(:,:,q),PhiQn(:,:,q),PhiFp(:,:,q),PhiFn(:,:,q),PsiQp(:,:,q),PsiFp(:,:,q),...
      PsiQn(:,:,q),PsiFn(:,:,q),lp(:,q),ln(:,q),s(:,q),kp(:,q),kn(:,q),~,PureN ]...
      = PolySolve_complex( w(q),Kc,Mc,Lc,nor,tol,lim,lim2);

end
ndof = ndofc;
l = Lc;
ci = v2struct (PhiQp,PhiQn,PhiFp,PhiFn,PsiQp,PsiFp,PsiQn,PsiFn,lp,ln,s,kp,kn,PureN,ndof,l);
clear PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn PureN
% [c] = sortDiff(ci,f);    
c = ci;

%% Should always use all propagating and nearfield waves, cant cut
nmodes_a = length(a.kp(:,1));
nmodes_b = length(b.kp(:,1));
nmodes_c = length(c.kp(:,1));

%     mode_limita = 5;
    mode_limita = nmodes_a;
    a.PhiQp = a.PhiQp(:,1:mode_limita,:);
    a.PhiQn = a.PhiQn(:,1:mode_limita,:);
    a.PsiQp = a.PsiQp(1:mode_limita,:,:);
    a.PsiQn = a.PsiQn(1:mode_limita,:,:);
    
    a.PhiFp = a.PhiFp(:,1:mode_limita,:);
    a.PhiFn = a.PhiFn(:,1:mode_limita,:);
    a.PsiFp = a.PsiFp(1:mode_limita,:,:);
    a.PsiFn = a.PsiFn(1:mode_limita,:,:);
    
%     mode_limitb = 3;
%     mode_limitb = length(knbi);
    mode_limitb = nmodes_b;
    b.PhiQp = b.PhiQp(:,1:mode_limitb,:);
    b.PhiQn = b.PhiQn(:,1:mode_limitb,:);
    
    b.PhiFp = b.PhiFp(:,1:mode_limitb,:);
    b.PhiFn = b.PhiFn(:,1:mode_limitb,:);
    
%     mode_limitc = 5;
     mode_limitc = nmodes_c;
    c.PhiQp = c.PhiQp(:,1:mode_limitc,:);
    c.PhiQn = c.PhiQn(:,1:mode_limitc,:);
    
    c.PhiFp = c.PhiFp(:,1:mode_limitc,:);
    c.PhiFn = c.PhiFn(:,1:mode_limitc,:);

%%

filename = ['Data/eigSolution_a' num2str(ndofa) 'b' num2str(ndofb) 'c' num2str(ndofc) ...
    'fi' num2str(fi) 'df' num2str(df) 'ff' num2str(ff)]
save(filename,'a', 'b', 'c','Ka','Kb','Kc','Ma','Mb','Mc','f');