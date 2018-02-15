% Find the reflection and transmission for three sections modeled in FE
% WFE Method
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

w = 2*pi*f;

%% Pre-allocate matrices

nmodes_a = length(a.kp(:,1));
nmodes_b = length(b.kp(:,1));
nmodes_c = length(c.kp(:,1));

lenf = length(f);
TRTr = zeros(1,lenf); TRTc = TRTr;

RPaa1 = zeros(nmodes_a,nmodes_a,lenf); 
TPba1 = zeros(nmodes_b,nmodes_a,lenf); 
RPbb1 = zeros(nmodes_b,nmodes_b,lenf); 
TPab1 = zeros(nmodes_a,nmodes_b,lenf); 

RPbb2 = zeros(nmodes_b,nmodes_b,lenf); 
TPcb2 = zeros(nmodes_c,nmodes_b,lenf); 
RPcc2 = zeros(nmodes_c,nmodes_c,lenf);
TPbc2 = zeros(nmodes_b,nmodes_c,lenf);

RPAA = RPaa1; TPCA = RPaa1;

for q=1:lenf 
     
    kPb = b.kp (1:nmodes_b,q);
    
    PsiQa_p0 = a.PsiQp(:,:,1);
    PsiFa_p0 = a.PsiFp(:,:,1);
    
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
    
end

TRTP = [RPAA; TPCA];

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
end

 %% Coefficient Plots

figure()
plot(f,abs(reshape(RPAA(1,1,:),[1 length(f)])),'b','LineWidth',2);

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|R|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

figure()
plot(f,abs(reshape(TPCA(1,1,:),[1 length(f)])),'b','LineWidth',2);

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','|T|','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])
  
%% Phase Plots

figure()
plot(f,phase(reshape(RPAA(1,1,:),[1 length(f)]))*180/pi,'b','LineWidth',2);

%  legend('Analytical Bending','WFE Bending', 'Analytical Longitudinal', 'WFE Longitudinal')
 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','\Phi_R','FontName','Times New Roman','FontSize',12)
set(gca,'fontsize',12,'FontName','Times New Roman');
%  axis([fi ff 0 2])

figure()
plot(f,phase(reshape(TPCA(1,1,:),[1 length(f)]))*180/pi,'b','LineWidth',2);

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

 set(get(gca,'XLabel'),'String','Frequency [Hz]','FontName','Times New Roman','FontSize',12)
 set(get(gca,'YLabel'),'String','Power coefficients','FontName','Times New Roman','FontSize',12)
 set(gca,'fontsize',12,'FontName','Times New Roman');
 
%% Save RT files

% fi = f(1);
% ff = f(end);
% df = mean(diff(f));
% 
% filename = ['RTW_' num2str(a.ndof) 'b' num2str(b.ndof) 'c' num2str(c.ndof) ...
%     'fi' num2str(fi) 'df' num2str(df) 'ff' num2str(ff)];
% save(filename, 'TRTP');