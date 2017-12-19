% Find the reflection and transmission for three sections modeled in FE
% Modular code:
% First create a eigsolution file!
% Breno Ebinuma Takiuti
% 12/14/2017

clear 
clc
% close all

%% Load Eig Solution

% load Data/eigSolution_a124b44c124fi100df1000ff250000
load Data/eigSolution_a16b8c16fi100df1000ff100000
% load Data/eigSolution_a16b8c16fi100df100ff100000
% load Data/eigSolution_a16b8c16fi100df100ff1000
% load Data/eigSolution_a16b8c16fi100df1ff1000

%% Geometric constants

L = 0.1;         % Use NeleB to calculate a length of B
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

md = 1;

for ii=1:NeleB
    Kabc( Inodab(:,:,ii), Inodab(:,:,ii)) = Kabc( Inodab(:,:,ii), Inodab(:,:,ii))+Kb;
    Mabc( Inodab(:,:,ii), Inodab(:,:,ii)) = Mabc( Inodab(:,:,ii), Inodab(:,:,ii))+Mb;
end
Kabc(Inodbc,Inodbc) = Kabc(Inodbc,Inodbc)+Kc;
Mabc(Inodbc,Inodbc) = Mabc(Inodbc,Inodbc)+Mc;

w = 2*pi*f;

for q=1:length(f)
 
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

    RA = eye(rp); RB = eye(rp);
    Z3 = zeros(rp);
    
    R = [RA Z3; Z3 RB];
    
%     Nind = [9,11,13,15];
    Nind = rp+1:2:rp*2;

    R(Nind,Nind) = -eye(length(Nind));


    D1 = (Kabc-w(q)^2*Mabc)*1;
    
   [rabc, cabc] = size(Kabc);
    
    InodE = [1:a.ndof/2,rabc-a.ndof/2+1:rabc];
    InodI = a.ndof/2+1:rabc-a.ndof/2;

%     %Method 1
    DEE = D1(InodE,InodE);
    DEI = D1(InodE,InodI);
    DIE = D1(InodI,InodE);
    DII = D1(InodI,InodI);
    
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

    RHAA(:,:,q) = Srt(1:cp,1:cp,q);
    THCA(:,:,q) = Srt(cp+1:cp+nModes,1:cp,q);
    
end

TRTH = [RHAA; THCA];

 %% Plots

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


%% Save RT files

filename = ['RTH_' num2str(ndofa) 'b' num2str(ndofb) 'c' num2str(ndofc) ...
    'fi' num2str(fi) 'df' num2str(df) 'ff' num2str(ff)];
save(filename, 'TRTH');