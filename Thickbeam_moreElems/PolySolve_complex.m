function [ PhiQa_p,PhiQa_n,PhiFa_p,PhiFa_n,PsiQa_p,PsiFa_p,PsiQa_n,PsiFa_n,llRpa,llRna,sa,kpa,kna,Phiq,pureN] = PolySolve_complex4( w,Ka,Ma,La,nor,tol,lim,lim2)
%% Using fan2016 as reference
%23/10/2017
% lim: used in pi/lim to calculate the tolerance to distinguish pure real
% and pure imaginary. bigger lim = tighter tolerance; smaller lim = looser
% tolerance
%
% lim2: used to distinguish between useful results and numerical modes.
% smaller is tighter
%
% tol: used to distinguish zeroes from very small values. Smaller is
% tighter

if nargin<7
    lim = 5;
    lim2 = 2000;
end

    [dofa,~] = size(Ka);
    dofa = dofa/2;

% PHIq and PHIf for the first section (a)
    % Define the dinamic matrix
    D = (Ka-w^2*Ma);
%     D = (Ka-w(q)^2*Ma);
    
    % Saparate it in quadrants
    DLL = (D(1:dofa,1:dofa));
    DLR = (D(1:dofa,dofa+1:dofa*2));
    DRL = (D(dofa+1:dofa*2,1:dofa));
    DRR = (D(dofa+1:dofa*2,dofa+1:dofa*2));
    
    Z=zeros(max(size(DRR)));
    I=eye(max(size(DRR)));
%     
%     A2 = DLR;
%     A1 = DLL+DRR;
%     A0 = DRL;

%     DD=[Z,I;-pinv(DLR)*DRL,-pinv(DLR)*(DLL+DRR)];
    
%     DD=[Z,I;-inv(DLR)*DRL,-inv(DLR)*(DLL+DRR)];   %changed this because
%     of the hybrid version. If it starts giving strange results. Restore.

    % Zhong
%     DD = [ -pinv(DLR)*DLL pinv(DLR); -DRL+DRR*pinv(DLR)*DLL -DRR*pinv(DLR)];
    %*Dont forget to change the Phif later
    
    % FAN and ICHCHOU 2016
    [~,N2] = size(DRR);
    N2 = N2^2;
    sig = norm(DRR,2)/N2;
    
    AA = [ Z sig*I; -DRL -DRR];
    BB = [ sig*I Z; DLL DLR];
    
    [VR,llR] = eig(AA,BB);
%     [VR,llR] = eigs(DD);
%     [VR,llR]=polyeig(A0,A1,A2);
    
    sa=1i.*log((diag(llR)));   %%llR=exp(-i*s*L) I find s
    sa=sa./La; %real positive value are for positive propagating
    % ee=10^(-4); %error
%     ee = tol;
    
    %% Normalize PhiQ!!!!
    
    [r,~] = size(VR);
    
    % Normalize with respect to the wave mode
%     if nor == 0
%         Phiq = VR;
%     else
%         Phiq = zeros(size(VR));
%         for ii = 1:r
%             Phiq(ii,:) = VR(ii,:)./(VR(nor,:));
%         end
%     end
%     llRd = diag(llR);
%     

    Phiq = zeros(size(VR));
    switch nor
        case {1,2}           
            for ii = 1:r
                Phiq(ii,:) = VR(ii,:)./(VR(nor,:));
            end
        case 3
            for ii = 1:r
                Phiq(ii,:) = VR(ii,:)./max(VR);
            end
        otherwise
            Phiq = VR;
    end
    llRd = diag(llR);
    %% S1: Separate only the pure real and pure imaginary wavenumbers
%     oc = unique_tol(sa,tol);
%     orderV = [];
    count = 1;
%     tol2 = 1e-4;
    tol2 = tol;
    
    % Generate a matrix of a combination of wavenumbers and their original
    % position on the wavenumbers vector: [k 1:1:end]_MODESx2
    s1 = [sa (1:length(sa))'];
    % Obs.: this is a trick to integrate the vector indexes into the
    % sorting process, so that any change in the vector's position is
    % "recorded" along the process. With this the second column of the
    % matrix can be used to "repeat the sorting" of the wavenumbers to the
    % wavemodes, without having to sort both along the process.
    
    s2 = s1;
    for ii = 1:length(sa)

        if (abs(real(sa(ii)*La))<pi/lim&&abs(imag(sa(ii)*La))<pi/lim)
            % Separate the numbers lower than pi over something
            if (abs(real(sa(ii)))>tol2 && abs(imag(sa(ii)))<tol2) || (abs(real(sa(ii)))<tol2 && abs(imag(sa(ii)))>tol2)
                % Separate pure real and pure imag and put then in the
                % beginning of the vector
                s2(count,:) = s1(ii,:);
                s2(ii,:) = s1(count,:);
                count = count+1;
            end
        end
    end
    nSort = count-1;
    
    % Sort the pure numbers:
    % First real then imaginary ('descend')
    warning('off','all')
%     s2
%     nSort
    [Pos, Neg] = pureRISort(s2(1:nSort,:),'descend');
    warning('on','all')
    
    % Remove the sorted number, retain only the complex (including errors)
    s5 = s2;
    s5(1:nSort,:) = [];
    
%     %% Take only the first smallest COMPLEX mode
%     
%     orderV = [];
%     s7 = [];
%     count = 1;
%     tol2 = tol;
%     for ii = 1:length(s5)
%         if (abs(real(s5(ii)))<400&&abs(imag(s5(ii)))<400)
%             s7(count,:) = s5(ii,:);
%             orderV(count)=ii;
%             count = count+1;
%         end
%     end
%     sa
%     s7
%     % Remove the complex small numbers, retain the rest
%     s5(orderV,:) = [];
%     
%     % Sort the complex numbers: first the real positive, if equal real,
%     % imaginary negative first
%     [Spos, Sneg] = complexSort2(s7,'ascend')
    
    %% Take only COMPLEX modes
    
    orderV = [];
    s8 = [];
    count = 1;
%     tol2 = 1e-4;
    tol2 = tol;
    for ii = 1:length(s5)
        if abs(real(s5(ii)))<lim2 && abs(real(s5(ii)))>tol2
            % this value of 400 is like pi/6: is an arbitrary value that
            % functions a separation between numerical errors and the valid results 
            s8(count,:) = s5(ii,:);
            orderV(count)=ii;
            count = count+1;
        end
    end
    s8;
    % Remove the complex small numbers, retain the rest
    s5(orderV,:) = [];

%     s8 = s5;

    % Sort the complex numbers: first the real positive, if equal real,
    % imaginary negative first
    warning('off','all')
    [Cpos, Cneg] = complexSort2(s8,'ascend',1e-2);
%     [Cpos, Cneg] = complexSort2(s8,'descend',1e-2);
    warning('on','all')
    size(Cpos);
    size(Cneg);
    %% Join back all the numbers:
    % First pure real, second the pure imaginary, third the complex and
    % last the errors
%     nRest = length(s5(:,1));
    
    % Sort the errors in the same way as the complex numbers, without
    % mixing with the other 3 types of numbers
%     [rp,rn] = complexSort(s5,'ascend');
    
    % Rebuild the wavenumbers vectors: separated in positive and negative
    % going waves
%     s6p = [Pos; Spos; Cpos; rp];
%     s6n = [Neg; Sneg; Cneg; rn];
%     s6p = [Pos];
%     s6n = [Neg];
    pureN = length(Pos);
    s6p = [Pos; Cpos];
    s6n = [Neg; Cneg];   
 
%     s6p = [Pos; Cpos; rp];
%     s6n = [Neg; Cneg; rn];

    % Separate the ordering vector
    index_p = real(s6p(:,2));
    index_n = real(s6n(:,2));
    
    [ppa,o] = size(s6p);
    [nn,o] = size(s6n);
    
    % Separate the wavenumbers vectors
    kpa = s6p(1:ppa,1);
    kna = s6n(1:nn,1);
    
    % Separate the lambda
    llRpa = llRd(index_p);
    llRna = llRd(index_n);
    
    % Calculate the "unsorted" wavemodes and left eigenvectors
    Phif = D*Phiq;
    Psiq = pinv(Phiq);
    Psif = pinv(Phif);
    
    Phiql1 = Phiq(1:dofa,:);
    Phiqa1 = Phiql1;
%     Phifa1 = Phif(1:dofa,:);
    Phiqr1 = Phiq(dofa+1:2*dofa,:);     %Zhong
    Phifa1 = DLL*Phiql1 + DLR*Phiqr1;
    Psiq1 = Psiq(:,1:dofa);
    Psif1 = Psif(:,1:dofa);
    
    % "Copy" the sorting used in the wavenumbers to the wavemodes and left
    % eigenvectors
    PhiQa_p = Phiqa1(:,index_p);
    PhiQa_n = Phiqa1(:,index_n);
    
    PsiQa_p = Psiq1(index_p,:);
    PsiQa_n = Psiq1(index_n,:);
    
    PhiFa_p = Phifa1(:,index_p);
    PhiFa_n = Phifa1(:,index_n);
    
    PsiFa_p = Psif1(index_p,:);
    PsiFa_n = Psif1(index_n,:);
    
%     a = struct('PhiQp',PhiQa_p,'PhiQn',PhiQa_n,'PhiFp',PhiFa_p,'PhiFn',PhiFa_n,'PsiQp',PsiQa_p,'PsiFp',PsiFa_p,'PsiQn',PsiQa_n,'PsiFn',PsiFa_n,'lp',llRpa,'ln',llRna,'s',sa,'kp',kpa,'kn',kna,'Phiq',Phiq,'PureN',pureN);
end

