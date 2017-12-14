function [PhiQ_pN, PhiQ_nN, PhiF_pN, PhiF_nN,kea] = WFE_reflection_beamRT_section(M,K,L,w)
%% Change in section (area)

dof = length(M)/2;  % total degrees of freedom in one side of the element

%% Eigenvalue problems
%% Polynomial eigenvalue problem
% known frequency vector (omega), find k

% Define the dinamic matrix
Da = (K-w^2*M);

% Saparate it in quadrants
DLL = (Da(1:dof,1:dof));
DLR = (Da(1:dof,dof+1:dof*2));
DRL = (Da(dof+1:dof*2,1:dof));
DRR = (Da(dof+1:dof*2,dof+1:dof*2));

% Define the polinomial eigenvalue problem
%   (lambda^2*DLR + lambda*(DLL+DRR) + DRL)qL = 0
%   (lambda^2*A2  + lambda*A1        + A0 )qL = 0
A2 = DLR;
A1 = DLL+DRR;
A0 = DRL;

% Solve the eigenvalue problem
% [Phia,lama] = polyeig(A0,A1,A2);
[Phia,lama,WR] = polyeig(A2,A1,A0);
lama;

kea = roundn((-log(lama)/L),-10);


%% Normalize PhiQ!!!!

[r,c] = size(Phia);

for ii = 1:r
    Phiqa(ii,:) = Phia(ii,:)./Phia(1,:);
end

% Phiqa = Phia;
Psia = WR';

%% Order the eigen solutions

ppa = 0; pna = 0; epa = 0; ena = 0;

for ii = 1:length(kea)
    Phifa(:,ii) = Da*[Phiqa(:,ii);lama(ii)*Phiqa(:,ii)];
    if isreal(lama(ii))==0 %propagating
        if imag(kea(ii))>0
            
            ppa = ppa+1;
            kpp(ppa) = kea(ii);  %positive propagating
            Phiqa_pp(:,ppa) = Phiqa(:,ii);
            Phifa_pp(:,ppa) = Phifa(1:dof,ii);
            
        elseif imag(kea(ii))<0
            
            pna = pna+1;
            kpn(pna) = kea(ii);  %negative propagating
            Phiqa_pn(:,pna) = Phiqa(:,ii);
            Phifa_pn(:,pna) = Phifa(1:dof,ii);
            
        end
    elseif isreal(lama(ii))==1 %evanescent
        if real(kea(ii))>0
            
            epa = epa+1;
            kep(epa) = kea(ii); %positive evanescent
            Phiqa_ep(:,epa) = Phiqa(:,ii);
            Phifa_ep(:,epa) = Phifa(1:dof,ii);
            
        elseif real(kea(ii))<0
            
            ena = ena+1;
            ken(ena) = kea(ii); %negative evanescent
            Phiqa_en(:,ena) = Phiqa(:,ii);
            Phifa_en(:,ena) = Phifa(1:dof,ii);
            
        end
    end
end

PhiQ_pN = [Phiqa_pp Phiqa_ep];

PhiQ_nN = [Phiqa_pn Phiqa_en];

PhiF_pN = [Phifa_pp Phifa_ep];

PhiF_nN = [Phifa_pn Phifa_en];

end

