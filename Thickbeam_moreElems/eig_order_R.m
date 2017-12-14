function [PhiQ_p,PhiQ_n,PhiF_p,PhiF_n,PsiQ_p,PsiQ_n,PsiF_p,PsiF_n,kppa1,kpna1,kepa1,kena1] = eig_order_R(dofa,D,s,Lx,Phiq,Phia1,Psiq1,Psif1,llR,ee)

% ee=10^(-1); %error

ppa = 0; pna = 0; epa = 0; ena = 0;
for ii = 1:length(s)
    Phifa1(:,ii) = D*[Phiq(:,ii)];
%     Phifa1(:,ii) = Phiq(dofa+1:2*dofa,ii);
    if abs(real(s(ii)*Lx))<pi/2&abs(imag(s(ii)*Lx))<pi/2
        
        if real(s(ii))>0+ee&(abs(imag(s(ii)))<0+ee)
            ppa = ppa+1;
            kppa1(ppa) = s(ii);  %positive propagating
            Phiq_pp1(:,ppa) = Phia1(:,ii);
            Phif_pp1(:,ppa) = Phifa1(1:dofa,ii);
            Psiq_pp1(ppa,:) = Psiq1(ii,:);
            Psif_pp1(ppa,:) = Psif1(ii,:);
            %             Phif_pp1(:,ppa) = Phifa1(:,ii);
            %Eigenvector to the right
            Phiq_ppR(:,ppa) = llR(ii,ii)*Phia1(:,ii);
            Phif_ppR(:,ppa) = Phifa1(dofa+1:2*dofa,ii);
            
        end
        if real(s(ii))<0-ee&(abs(imag(s(ii)))<0+ee)
            pna = pna+1;
            kpna1(pna) = s(ii);  %negative propagating
            Phiq_pn1(:,pna) = Phia1(:,ii);
            Phif_pn1(:,pna) = Phifa1(1:dofa,ii);
            Psiq_pn1(pna,:) = Psiq1(ii,:);
            Psif_pn1(pna,:) = Psif1(ii,:);
%             Phif_pn1(:,pna) = Phifa1(:,ii);
            %Right 
            Phiq_pnR(:,pna) = llR(ii,ii)*Phia1(:,ii);
            Phif_pnR(:,pna) = Phifa1(dofa+1:2*dofa,ii);
        end
        if imag(s(ii))<0-ee&(abs(real(s(ii)))<0+ee)
            epa = epa+1;
            kepa1(epa) = s(ii); %positive evanescent
            Phiq_ep1(:,epa) = Phia1(:,ii);
            Phif_ep1(:,epa) = Phifa1(1:dofa,ii);
            Psiq_ep1(epa,:) = Psiq1(ii,:);
            Psif_ep1(epa,:) = Psif1(ii,:);
%             Phif_ep1(:,epa) = Phifa1(:,ii);
            %R
            Phiq_epR(:,epa) = llR(ii,ii)*Phia1(:,ii);
            Phif_epR(:,epa) = Phifa1(dofa+1:2*dofa,ii);
        end
        if imag(s(ii))>0+ee&(abs(real(s(ii)))<0+ee)
            ena = ena+1;
            kena1(ena) = s(ii); %negative evanescent
            Phiq_en1(:,ena) = Phia1(:,ii);
            Phif_en1(:,ena) = Phifa1(1:dofa,ii);
            Psiq_en1(ena,:) = Psiq1(ii,:);
            Psif_en1(ena,:) = Psif1(ii,:);
%             Phif_en1(:,ena) = Phifa1(:,ii);
            %R
            Phiq_enR(:,ena) = llR(ii,ii)*Phia1(:,ii);
            Phif_enR(:,ena) = Phifa1(dofa+1:2*dofa,ii);
        end
    end
end

%     kp = [kppa1 kepa1];
%     kn = [kpna1 kena1];
kp = [kppa1 kepa1];
kn = [kpna1 kena1];

% PhiQ_p = [Phiq_pp1 Phiq_ep1];
PhiQ_p = [Phiq_ppR Phiq_epR];
PsiQ_p = [Psiq_pp1; Psiq_ep1];
% PhiQ_p = [Phiq_pp1];
% PhiQ_p = fliplr(PhiQ_p);

% PhiQ_n = [Phiq_pn1 Phiq_en1];
PhiQ_n = [Phiq_pnR Phiq_enR];
PsiQ_n = [Psiq_pn1; Psiq_en1];
% PhiQ_n = [Phiq_pn1];
% PhiQ_n = fliplr(PhiQ_n);

% PhiF_p = [Phif_pp1 Phif_ep1];
PsiF_p = [Psif_pp1; Psif_ep1];
PhiF_p = [Phif_ppR Phif_epR];
% PhiF_p = [Phif_pp1];
% PhiF_p = fliplr(PhiF_p);

% PhiF_n = [Phif_pn1 Phif_en1];
PsiF_n = [Psif_pn1; Psif_en1];
PhiF_n = [Phif_pnR Phif_enR];
% PhiF_n = [Phif_pn1];
% PhiF_n = fliplr(PhiF_n);