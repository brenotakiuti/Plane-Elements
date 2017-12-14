function [PhiQ_pN, PhiQ_nN, PhiF_pN, PhiF_nN,PsiQ_pN,PsiQ_nN,PsiF_pN,PsiF_nN,kea] = WFE_reflection_beamRT_section2(Ma,Ka,Lx,w)
%% Change in section (area)

dofa = length(Ma)/2;

% PHIq and PHIf for the first section (a)
% Define the dinamic matrix
D = (Ka-w^2*Ma);

% Saparate it in quadrants
DLL = (D(1:dofa,1:dofa));
DLR = (D(1:dofa,dofa+1:dofa*2));
DRL = (D(dofa+1:dofa*2,1:dofa));
DRR = (D(dofa+1:dofa*2,dofa+1:dofa*2));

% Z=zeros(max(size(DRR)));
%     I=eye(max(size(DRR)));
%     DD=[Z,I;-inv(DLR)*DRL,-inv(DLR)*(DLL+DRR)];
%
%% Other formulation, Renno (2013)
DD = [ -DLR^-1*DLL DLR^-1; -DRL+DRR*DLR^-1*DLL -DRR*DLR^-1];    % Renno(2013)
[VR,llR,WR]=eig(DD);

% Phia=VR(1:dofa,:);

%VR(9:16,1)=VR(1:8,1)*llR(1,1)  they represent the eigenvectors for disp at x=L and x=R
%         [VL,llL]=eig(DD'); %conj transpose of the left eigenvector
%         VL=VL';
%         [ll,o]=sort(diag(llR));
%            for g=1:max(size(llR))
%                VRs(:,g)=VR(:,o(g));
%            end
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=i.*log((diag(llR)));   %%llR=exp(-i*s*L) I find s
s=s./Lx; %real positive value are for positive propagating
%k_eli(:,jj)=s;
ee=10^(-4); %error

%% Normalize PhiQ!!!!

[r,c] = size(VR);

% Phiq = Phi;
for ii = 1:r
    Phiq(ii,:) = VR(ii,:)./VR(1,:);
end

pos = [8 11 12 15]; %500Hz
% pos = [11 12 13 14]; %1000Hz
% Phiq(1:2:dofa*2,pos) = 0;
% Psiq = inv(Phiq);

Psiq = WR';
%%
Phia1=Phiq(1:dofa,:);
Psiq1 = Psiq(:,1:dofa);
Psif1 = Psiq(:,dofa+1:2*dofa);

%%
ppa = 0; pna = 0; epa = 0; ena = 0;
for ii = 1:length(s)
%     Phifa1(:,ii) = D*[Phiq(:,ii)];
    Phifa1(:,ii) = Phiq(dofa+1:2*dofa,ii);
    if abs(real(s(ii)*Lx))<pi/10&abs(imag(s(ii)*Lx))<pi/10
        
        if real(s(ii))>0+ee&(abs(imag(s(ii)))<0+ee)
            ppa = ppa+1;
            kppa1(ppa) = s(ii);  %positive propagating
            Phiq_pp1(:,ppa) = Phia1(:,ii);
%             Phif_pp1(:,ppa) = Phifa1(1:dofa,ii);
            Phif_pp1(:,ppa) = Phifa1(:,ii);
            Psiq_pp1(ppa,:) = Psiq1(ii,:);
            Psif_pp1(ppa,:) = Psif1(ii,:);
        end
        if real(s(ii))<0-ee&(abs(imag(s(ii)))<0+ee)
            pna = pna+1;
            kpna1(pna) = s(ii);  %negative propagating
            Phiq_pn1(:,pna) = Phia1(:,ii);
%             Phif_pn1(:,pna) = Phifa1(1:dofa,ii);
            Phif_pn1(:,pna) = Phifa1(:,ii);
            Psiq_pn1(pna,:) = Psiq1(ii,:);
            Psif_pn1(pna,:) = Psif1(ii,:);
        end
        if imag(s(ii))<0-ee&(abs(real(s(ii)))<0+ee)
            epa = epa+1;
            kepa1(epa) = s(ii); %positive evanescent
            Phiq_ep1(:,epa) = Phia1(:,ii);
%             Phif_ep1(:,epa) = Phifa1(1:dofa,ii);
            Phif_ep1(:,epa) = Phifa1(:,ii);
            Psiq_ep1(epa,:) = Psiq1(ii,:);
            Psif_ep1(epa,:) = Psif1(ii,:);
        end
        if imag(s(ii))>0+ee&(abs(real(s(ii)))<0+ee)
            ena = ena+1;
            kena1(ena) = s(ii); %negative evanescent
            Phiq_en1(:,ena) = Phia1(:,ii);
%             Phif_en1(:,ena) = Phifa1(1:dofa,ii);
            Phif_en1(:,ena) = Phifa1(:,ii);
            Psiq_en1(ena,:) = Psiq1(ii,:);
            Psif_en1(ena,:) = Psif1(ii,:);
        end
    end
end

%     kp = [kppa1 kepa1];
%     kn = [kpna1 kena1];
kp = [kppa1];
kn = [kpna1];

PhiQ_pN = [Phiq_pp1 Phiq_ep1];
PsiQ_pN = [Psiq_pp1; Psiq_ep1];
% PhiQ_p = [Phiq_pp1];
% PhiQ_p = fliplr(PhiQ_p);

PhiQ_nN = [Phiq_pn1 Phiq_en1];
PsiQ_nN = [Psiq_pn1; Psiq_en1];
% PhiQ_n = [Phiq_pn1];
% PhiQ_n = fliplr(PhiQ_n);

PhiF_pN = [Phif_pp1 Phif_ep1];
PsiF_pN = [Psif_pp1; Psif_ep1];
% PhiF_p = [Phif_pp1];
% PhiF_p = fliplr(PhiF_p);

PhiF_nN = [Phif_pn1 Phif_en1];
PsiF_nN = [Psif_pn1; Psif_en1];
% PhiF_n = [Phif_pn1];
% PhiF_n = fliplr(PhiF_n);

end

