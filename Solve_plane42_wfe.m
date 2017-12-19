% Find the reflection and transmission for three sections modeled in FE
% Modular code:
% First create a eigsolution file!
% Breno Ebinuma Takiuti
% 12/14/2017

clear 
clc
% close all

%% Load Eig Solution

load Data/eigSolution_a124b44c124fi100df1000ff100000
% load Data/eigSolution_a52b20c52fi100df1000ff100000
% load Data/eigSolution_a124b44c124fi100df1000ff250000
% load Data/eigSolution_a16b8c16fi100df1000ff100000
% load Data/eigSolution_a16b8c16fi100df100ff100000
% load Data/eigSolution_a16b8c16fi100df100ff1000
% load Data/eigSolution_a16b8c16fi100df1ff1000

%% Geometric constants

L = b.l;
% L = 0.1;         % Use NeleB to calculate a length of B
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

for ii=1:NeleB
    Inodab(:,:,ii) = [round(a.ndof/2*ii)+Inoda, a.ndof+(a.ndof/2*(ii-1))+Inoda];
end
Inodbc = [a.ndof+1+(a.ndof/2*(NeleB-1)):a.ndof+(a.ndof/2*(NeleB+1))];

Ea = eye(a.ndof/2);
Eb(Inoda,Inodb) = I;
Ec = Ea;
Ca(Inoda,Inoda) = I;
Cb = Eb;
Cc = Ca;

w = 2*pi*f;

for q=1:length(f)
    
    %% Numeric
     
    nmodes_b = length(b.kp(:,q));
    nmodes_a = length(a.kp(:,q));
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

 %% Cutoff plots
% Use the plotCutOffs only for the non zero coefficients, the others plot
% them normally
% 
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
 
%% Save RT files

fi = f(1);
ff = f(end);
df = mean(diff(f));

filename = ['Data/RTW_' num2str(a.ndof) 'b' num2str(b.ndof) 'c' num2str(c.ndof) ...
    'fi' num2str(fi) 'df' num2str(df) 'ff' num2str(ff)];
save(filename, 'TRTP');