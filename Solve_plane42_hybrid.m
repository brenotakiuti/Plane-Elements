% Find the reflection and transmission for three sections modeled in FE
% HYBRID method, by RENNO and MITROU
% Modular code:
% First create a eigsolution file!
% Breno Ebinuma Takiuti
% 12/14/2017
tic

clear 
clc
% close all

%% Load Eig Solution

% load Data/eigSolution_a52b20c52fi100df100ff130000
% load Data/eigSolution_a124b44c124fi100df1000ff250000
% load Data/eigSolution_a16b8c16fi100df100ff130000
% load Data/eigSolution_a16b8c16fi100df1000ff135000
load Data/eigSolution_a16b8c16fi100df1000ff100000
% load Data/eigSolution_a16b8c16fi100df100ff100000
% load Data/eigSolution_a16b8c16fi100df100ff1000
% load Data/eigSolution_a16b8c16fi100df1ff1000

%% Geometric constants

n = 10;
% L = 0.1;         % Use NeleB to calculate a length of B
L = n*b.l;
NeleB = round(L/b.l);            % Number of elements in B

nModes = 2;

%% Joint FE model

Kabc = zeros(a.ndof+(a.ndof/2*(NeleB+1)));
Mabc = Kabc;

Kabc(1:a.ndof,1:a.ndof) = Ka;
Mabc(1:a.ndof,1:a.ndof) = Ma;

%% Boundary conditions when using the left eigenvectors

Eb = zeros(a.ndof/2,b.ndof/2);
Ca = zeros(a.ndof/2);

% Describe the boundary (connection a-b)
I = eye(b.ndof/2);
Inoda = (a.ndof/2-b.ndof/2)/2+1:(a.ndof/2-b.ndof/2)/2+b.ndof/2; %SM
%     Inoda = (ar/2-44/2)/2+1:(ar/2-44/2)/2+br/2; %Quasi SM
%     Inoda = 1:br/2;                     %NS
%     Inoda = (ar/2-br/2)+1:ar/2;
% Inoda = (ar-(ra/2)-br/2)/2+1:(ar-(ra/2)-br/2)/2+br/2; %SM2
Inodb = 1:b.ndof/2;

Inodab = zeros(1,length(Inoda)*2,NeleB);
for ii=1:NeleB
    Inodab(:,:,ii) = [round(a.ndof/2*ii)+Inoda, a.ndof+(a.ndof/2*(ii-1))+Inoda];
end
Inodbc = a.ndof+1+(a.ndof/2*(NeleB-1)):a.ndof+(a.ndof/2*(NeleB+1));

Ea = eye(a.ndof/2);
Eb(Inoda,Inodb) = I;
Ec = Ea;
Ca(Inoda,Inoda) = I;
Cb = Eb;
Cc = Ca;

for ii=1:NeleB
    Kabc( Inodab(:,:,ii), Inodab(:,:,ii)) = Kabc( Inodab(:,:,ii), Inodab(:,:,ii))+Kb;
    Mabc( Inodab(:,:,ii), Inodab(:,:,ii)) = Mabc( Inodab(:,:,ii), Inodab(:,:,ii))+Mb;
end
Kabc(Inodbc,Inodbc) = Kabc(Inodbc,Inodbc)+Kc;
Mabc(Inodbc,Inodbc) = Mabc(Inodbc,Inodbc)+Mc;

w = 2*pi*f;

%% Pre-allocate matrices

[ndof_a,nmodes_a,~] = size(a.PhiQp);
[ndof_b,nmodes_b,~] = size(b.PhiQp);
Z1 = zeros(size(a.PhiQp(:,:,1)));
Z2 = Z1';
I = eye(size(a.PhiQp(:,:,1)));
lenf = length(f);

Srt = zeros(nmodes_a*2,nmodes_a*2,lenf);
RHAA = zeros(nmodes_a,nmodes_a,lenf);
THCA = RHAA;

RA = eye(ndof_a); 
RB = eye(ndof_a);
Z3 = zeros(ndof_a);
R = [RA Z3; Z3 RB];
[rabc, cabc] = size(Kabc);
    
InodE = [1:a.ndof/2,rabc-a.ndof/2+1:rabc];  % External DOFs
InodI = a.ndof/2+1:rabc-a.ndof/2;           % Internal DOFs

for q=1:lenf
 
    %% Hybrid method by Renno

    PhiQ_p = [a.PhiQp(:,:,q) Z1; Z1 c.PhiQp(:,:,q)];
    PhiQ_n = [a.PhiQn(:,:,q) Z1; Z1 c.PhiQn(:,:,q)];
    PhiF_p = [a.PhiFp(:,:,q) Z1; Z1 c.PhiFp(:,:,q)];
    PhiF_n = [a.PhiFn(:,:,q) Z1; Z1 c.PhiFn(:,:,q)];
    

    PsiQ_n = [a.PsiQn(:,:,q) Z2; Z2 c.PsiQn(:,:,q)];

    Nind = ndof_a+1:2:ndof_a*2;

    % Generate a rotation matrix
    R(Nind,Nind) = -eye(length(Nind));

    D0 = (Kabc-w(q)^2*Mabc);

%     %Method 1
    DEE = D0(InodE,InodE);
    DEI = D0(InodE,InodI);
    DIE = D0(InodI,InodE);
    DII = D0(InodI,InodI);
    
    Dj = (DEE-DEI*pinv(DII)*DIE);    %Method 1

    %Method 2
%     ni = length(InodE);
%     DTLL = D1(InodE(1:ni/2),InodE(1:ni/2));
%     DTLR = D1(InodE(1:ni/2),InodE(ni/2+1:ni));
%     DTRL = D1(InodE(ni/2+1:ni),InodE(1:ni/2));
%     DTRR = D1(InodE(ni/2+1:ni),InodE(ni/2+1:ni));
%     DTLO = D1(InodE(1:ni/2),InodI);
%     DTRO = D1(InodE(ni/2+1:ni),InodI);
%     DTOL = D1(InodI,InodE(1:ni/2));
%     DTOR = D1(InodI,InodE(ni/2+1:ni));
%     DTOO = D1(InodI,InodI);
%     
%     DLL = DTLL - DTLO*pinv(DTOO)*DTOL;
%     DLR = DTLR - DTLO*pinv(DTOO)*DTOR;
%     DRL = DTRL - DTRO*pinv(DTOO)*DTOL;
%     DRR = DTRR - DTRO*pinv(DTOO)*DTOR;
% %     
%     Dj = [DLL DLR; DRL DRR];    %Method 2

    %Method 3
%     Dj = D0;
        
    Srt(:,:,q) = pinv(PsiQ_n*(-Dj*R*PhiQ_n+R*PhiF_n))*PsiQ_n*(Dj*R*PhiQ_p-R*PhiF_p);
%     Srt(:,:,q) = -pinv(Dj*R*PhiQ_n-R*PhiF_n)*(Dj*R*PhiQ_p-R*PhiF_p);

    RHAA(:,:,q) = Srt(1:nmodes_a,1:nmodes_a,q);
    THCA(:,:,q) = Srt(nmodes_a+1:end,1:nmodes_a,q);
    
end

TRTH = [RHAA; THCA];

toc

%% Calculate Power coefficients

PrPP2 = zeros(1,lenf);
PrPL2 = PrPP2; PrPN2 = PrPP2; PtPP2 = PrPP2; PtPL2 = PrPP2; PtPN2 = PrPP2;
PrLP2 = PrPP2; PrLL2 = PrPP2; PrLN2 = PrPP2; PtLP2 = PrPP2; PtLL2 = PrPP2; PtLN2 = PrPP2; 

for q=1:lenf 
    % Power Matrix
    % MITROU (2015)
    Pa2 = (1i*w(q)/2)*[a.PhiQp(:,:,q)'*a.PhiFp(:,:,q) a.PhiQp(:,:,q)'*a.PhiFn(:,:,q);
        a.PhiQn(:,:,q)'*a.PhiFp(:,:,q) a.PhiQn(:,:,q)'*a.PhiFn(:,:,q)]-...
        [a.PhiFp(:,:,q)'*a.PhiQp(:,:,q) a.PhiFp(:,:,q)'*a.PhiQn(:,:,q);
        a.PhiFn(:,:,q)'*a.PhiQp(:,:,q) a.PhiFn(:,:,q)'*a.PhiQn(:,:,q)];
    Pc2 = (1i*w(q)/2)*[c.PhiQp(:,:,q)'*c.PhiFp(:,:,q) c.PhiQp(:,:,q)'*c.PhiFn(:,:,q);
        c.PhiQn(:,:,q)'*c.PhiFp(:,:,q) c.PhiQn(:,:,q)'*c.PhiFn(:,:,q)]-...
        [c.PhiFp(:,:,q)'*c.PhiQp(:,:,q) c.PhiFp(:,:,q)'*c.PhiQn(:,:,q);
        c.PhiFn(:,:,q)'*c.PhiQp(:,:,q) c.PhiFn(:,:,q)'*c.PhiQn(:,:,q)];
    
    % Power Coefficients
    PrPP2(q) = abs(RHAA(1,1,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(1,1));
    PrPL2(q) = abs(RHAA(2,1,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(1,1));
    PrPN2(q) = abs(RHAA(2,1,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(1,1));
%     PrPC12(q) = abs(RPAA(4,1,q))^2*(Pa2(nmodes_a+4,nmodes_a+4)/Pa2(1,1));
%     PrPC22(q) = abs(RPAA(5,1,q))^2*(Pa2(nmodes_a+5,nmodes_a+5)/Pa2(1,1));
    PtPP2(q) = abs(THCA(1,1,q))^2*(Pc2(1,1)/Pa2(1,1));
    PtPL2(q) = abs(THCA(2,1,q))^2*(Pc2(2,2)/Pa2(1,1));
    PtPN2(q) = abs(THCA(3,1,q))^2*(Pc2(3,3)/Pa2(1,1));
%     PtPC12(q) = abs(TPCA(4,1,q))^2*(Pc2(4,4)/Pa2(1,1));
%     PtPC22(q) = abs(TPCA(5,1,q))^2*(Pc2(5,5)/Pa2(1,1));
   
    PrLP2(q) = abs(RHAA(1,2,q))^2*(Pa2(nmodes_a+1,nmodes_a+1)/Pa2(2,2));
    PrLL2(q) = abs(RHAA(2,2,q))^2*(Pa2(nmodes_a+2,nmodes_a+2)/Pa2(2,2));
    PrLN2(q) = abs(RHAA(3,2,q))^2*(Pa2(nmodes_a+3,nmodes_a+3)/Pa2(2,2));
%     PrLC12(q) = abs(RPAA(4,2,q))^2*(Pa2(nmodes_a+4,nmodes_a+4)/Pa2(2,2));
%     PrLC22(q) = abs(RPAA(5,2,q))^2*(Pa2(nmodes_a+5,nmodes_a+5)/Pa2(2,2));
    PtLP2(q) = abs(THCA(1,2,q))^2*(Pc2(1,1)/Pa2(2,2));
    PtLL2(q) = abs(THCA(2,2,q))^2*(Pc2(2,2)/Pa2(2,2));
    PtLN2(q) = abs(THCA(3,2,q))^2*(Pc2(3,3)/Pa2(2,2));
%     PtLC12(q) = abs(TPCA(4,2,q))^2*(Pc2(4,4)/Pa2(2,2));
%     PtLC22(q) = abs(TPCA(5,2,q))^2*(Pc2(5,5)/Pa2(2,2));
end

 %% Coefficient Plots

% 
figure()
plot(f,abs(reshape(RHAA(1,1,:),[1 length(f)])),'r--')

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

figure()
plot(f,abs(reshape(THCA(1,1,:),[1 length(f)])),'r--')

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

 %% Phase Plots

% 
figure()
plot(f,phase(reshape(RHAA(1,1,:),[1 length(f)]))*180/pi,'r--')

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','\Phi_R','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

figure()
plot(f,phase(reshape(THCA(1,1,:),[1 length(f)]))*180/pi,'r--')

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','\Phi_T','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

%% Power Plots
 figure()
 plot(f,abs(PrPP2),'b','LineWidth',1)
 hold on
%  plot(f,abs(PrPL2),'r','LineWidth',1)
 plot(f,abs(PrPN2),'g','LineWidth',1)
%  plot(f,abs(PrPC12),'m','LineWidth',1)
%  plot(f,abs(PrPC22),'c','LineWidth',1)
 plot(f,abs(PtPP2),'b--','LineWidth',1)
%  plot(f,abs(PtPL2),'r*','LineWidth',1)
 plot(f,abs(PtPN2),'g--','LineWidth',1)
%  plot(f,abs(PtPC12),'m*','LineWidth',1)
%  plot(f,abs(PtPC22),'c*','LineWidth',1)

% for ii = 1:length(Cutoff_a)
%     plot([Cutoff_ai(ii) Cutoff_ai(ii)],[-800 800],'k:','LineWidth',0.5)
% end
% for ii = 1:length(Cutoff_b)
%     plot([Cutoff_bi(ii) Cutoff_bi(ii)],[-800 800],'k:','LineWidth',0.5)
% end

%  legend('Propagating bending reflected','Longitudinal reflected','Nearfield reflected','Attenuated 1 reflected','Attenuated 2 reflected','Propagating bending transmitted','Longitudinal transmitted','Nearfield transmitted','Attenuated 1 transmitted','Attenuated 2 transmitted')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','Power coefficients','FontName','Times New Roman','FontSize',12)
 set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi f(q) -0.01 2])

%% Save RT files
% 
% fi = f(1);
% ff = f(end);
% df = mean(diff(f));
% 
% filename = ['RTH_' num2str(a.ndof) 'b' num2str(b.ndof) 'c' num2str(c.ndof) ...
%     'fi' num2str(fi) 'df' num2str(df) 'ff' num2str(ff)];
% save(filename, 'TRTH');