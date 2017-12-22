% Find the reflection and transmission for three sections modeled in FE
% Modular code:
% First create a eigsolution file!
% Breno Ebinuma Takiuti
% 12/14/2017

clear 
clc
% close all

%% Load Eig Solution

load Data/eigSolution_a52b20c52fi100df100ff130000
% load Data/eigSolution_a124b44c124fi100df1000ff250000
% load Data/eigSolution_a16b8c16fi100df100ff130000
% load Data/eigSolution_a16b8c16fi100df1000ff135000
% load Data/eigSolution_a16b8c16fi100df1000ff100000
% load Data/eigSolution_a16b8c16fi100df100ff100000
% load Data/eigSolution_a16b8c16fi100df100ff1000
% load Data/eigSolution_a16b8c16fi100df1ff1000

%% Material's Constants
% Material: Steel

rho=7800;     %mass per unit valume
E=2.06e11;  %Young's modulus

%% Geometric constants

tb = a.l;              % base of the cross-section (tickness) (m)
ha = 7.5e-3;              % height of the cross-section (width) (m)
hb = 2.56e-3;

% ha = 18e-3;
% hb = 6e-3;

hc = ha;

La = a.l;
Lb = b.l;               % length of the element (x direction) (m) (3 elem)
Lc = c.l;

%% Discontinuity size
n = 100;
% L = 0.1;         % Use NeleB to calculate a length of B
L = n*b.l;

%%
Sa = tb*(ha);
Sb = tb*hb;
Sc = tb*hc;
Ia = tb*(ha)^3/12;
Ib = tb*(hb)^3/12;
Ic = tb*(hc)^3/12;
beta_ab = Sb/Sa;
beta_bc = Sc/Sb;

nModes = 2;

w = 2*pi*f;

for q=1:length(f)
    %% Analytical Solution
    
    AA = eye(2);
    BB = AA;

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
     
        
end

% TRTA = [RBTAA; TBTCA];

%% Coefficient Plots

figure()
plot(f,abs(reshape(RBTAA(1,1,:),[1 length(f)])),'b:')

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

figure()
plot(f,abs(reshape(TBTCA(1,1,:),[1 length(f)])),'b:')

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

%% Phase Plots

figure()
plot(f,phase(reshape(RBTAA(1,1,:),[1 length(f)])),'b:')

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

figure()
plot(f,phase(reshape(TBTCA(1,1,:),[1 length(f)])),'b:')

%  plot(f,abs(R_WFE3(3,:)),'m-.','LineWidth',3)
%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

%% Save RT files

% fi = f(1);
% ff = f(end);
% df = mean(diff(f));
% filename = ['RTA_' num2str(ndofa) 'b' num2str(ndofb) 'c' num2str(ndofc) ...
%     'fi' num2str(fi) 'df' num2str(df) 'ff' num2str(ff)];
% save(filename, 'TRTA');