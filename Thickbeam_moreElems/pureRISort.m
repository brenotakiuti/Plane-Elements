function [Pos,Neg]=pureRISort(X,MODE,tol)
% Works symilarly to the command sort but sorts complex numbers considering
% the real part as always bigger than the imaginary. Ex.:
%               1+1i>0+2i
%               1+1i>0+1.4142i
%               1+0i>0+1000i
% If the 'descent'  mode is chosen, the numbers with bigger real parts will
% always be in the first half of the sorted vector.
% Only works for column vectors

if nargin<4
    if nargin<3
        ee = 1e-3;
        if nargin<2
            MODE = 'ascend';
        end
    else
        ee = tol;
    end
else
    error('number of arguments is invalid')
end

n = length(X);
Y = zeros(n,1);
J = Y;

% XI = [X I];

count_R = 1;
count_I = 1;

R_part = [0 0];
I_part = [0 0];

for ii=1:n
    if abs(real(X(1,1)))>ee && abs(imag(X(1,1)))<ee
        % Is REAL
        R_part(count_R,:) = X(1,:);
        count_R = count_R+1;
        X(1,:) = [];
    elseif abs(imag(X(1,1)))>ee && abs(real(X(1,1)))<ee
        % Is IMAGINARY
        I_part(count_I,:) = X(1,:);
        count_I = count_I+1;
        X(1,:) = [];
    else
        error(['Complex value found, in position [' num2str(ii) '] either eliminate the complex values or loosen the tolerance']);
    end
end

[nR,~] = size(R_part);
[nI,~] = size(I_part);

if strcmp(MODE,'ascend') || strcmp(MODE,'Ascend')
    
    %Sort real
    for ii=1:nR-1
        for jj=1:nR-1
            if real(R_part(jj+1))<real(R_part(jj))
                [R_part]=swap(R_part,jj);
            end
        end
    end
    
    %Sort imag
    for ii=1:nI-1
        for jj=1:nI-1
            if imag(I_part(jj+1))<imag(I_part(jj))
                [I_part]=swap(I_part,jj);
            end
        end
    end
    
elseif strcmp(MODE,'descend') || strcmp(MODE,'Descend')
    
    %Sort real
    for ii=1:nR-1
        for jj=1:nR-1
            if real(R_part(jj+1))>real(R_part(jj))
                [R_part]=swap(R_part,jj);
            end
        end
    end
    
    %Sort imag
    for ii=1:nI-1
        for jj=1:nI-1
            if imag(I_part(jj+1))>imag(I_part(jj))
                [I_part]=swap(I_part,jj);
            end
        end
    end
    
else
    error('Invalid value for MODE, please choose "ascend" or "descend"');
end

    R_part(nR/2+1:end,:) = flipud(R_part(nR/2+1:end,:));
    I_part(nI/2+1:end,:) = flipud(I_part(nI/2+1:end,:));
    
    Rp = R_part(1:nR/2,:);
    Rn = R_part(nR/2+1:end,:);
    Ip = I_part(1:nI/2,:);
    In = I_part(nI/2+1:end,:);
    
    Pos = [Rp; In] ;
    Neg = [Rn; Ip] ;
end

function [X]=swap(X,ii)
swp = X(ii,:);
X(ii,:)=X(ii+1,:);
X(ii+1,:)=swp;
% swp = I(ii);
% I(ii)=I(ii+1);
% I(ii+1)=swp;
end