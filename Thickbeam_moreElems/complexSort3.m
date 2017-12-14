function [Pos,Neg]=complexSort3(X,MODE,tol)
% Works symilarly to the command sort but sorts complex numbers considering
% the real part as always bigger than the imaginary. Ex.:
%               1+1i>0+2i
%               1+1i>0+1.4142i
%               1+0i>0+1000i
% If the 'descent'  mode is chosen, the numbers with bigger real parts will
% always be in the first half of the sorted vector.
% Only works for column vectors
Pos = [];
Neg = [];

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

% [X(:,1),in] = sort(X(:,1),'descend');
% X(:,2:end) = X(in,2:end);

[n,~] = size(X);
% ee = 1e-3;

Y = zeros(size(X));

if strcmp(MODE,'descend') || strcmp(MODE,'Descend')
    for ii=1:n-1
        for jj=1:n-1
            if real(X(jj,1)) < real(X(jj+1,1))
                X=swap(X,jj);
            end
            if abs(real(X(jj,1)) - real(X(jj+1,1)))<ee
                if imag(X(jj+1,1)) < imag(X(jj,1))
                    X=swap(X,jj);
                end
            end
        end
    end
    Pos = (X(1:n/2,:));
    Negi = (X(n/2+1:end,:));
    
    [nNeg,~] = size(Negi);
    
    for ii = 0:2:nNeg-1
        ii;
        Neg(ii+1:ii+2,:) = Negi(nNeg-ii-1:nNeg-ii,:);
    end
    
elseif strcmp(MODE,'ascend') || strcmp(MODE,'Ascend')
    countp = 0;
    countn = 0;
    for jj=1:n
        if (real(X(jj,1)) > 0 && imag(X(jj,1))>0) || (real(X(jj,1)) < 0 && imag(X(jj,1))<0)
            countp = countp+1;
            Px(countp,:) = X(jj,:);
            
        else
            countn = countn+1;
            Nx(countn,:) = X(jj,:);
            
        end
    end
%     X
%     Px
%     countn
%     countp

    [Px1]=sortiasc(Px,countp,ee);
    [Nx1]=sortiasc(Nx,countn,ee);
    
    Pos = Px1;
    Neg = Nx1;
    
%     [nNeg,~] = size(Negi);
    
%     for ii = 0:2:nNeg-1
%         ii;
%         Neg(ii+1:ii+2,:) = Negi(nNeg-ii-1:nNeg-ii,:);
%     end
    
else
    error('Invalid value for MODE, please choose "ascend" or "descend"');
end

X;
end

function [X]=swap(X,ii)
swp = X(ii,:);
X(ii,:)=X(ii+1,:);
X(ii+1,:)=swp;
% swp = I(ii);
% I(ii)=I(ii+1);
% I(ii+1)=swp;
end

function [Y]=sortiasc(X,countn,ee)
Y = zeros(size(X));
for ii=1:countn-1
    for jj=1:countn-1
        if real(X(jj,1)) > real(X(jj+1,1))
            X=swap(X,jj);
        end
        if abs(real(X(jj,1)) - real(X(jj+1,1)))<ee
            if imag(X(jj+1,1)) < imag(X(jj,1))
                X=swap(X,jj);
            end
        end
    end
end
Y(1:2:end,:) = (X(countn/2+1:end,:));
Y(2:2:end,:) = flipud(X(1:countn/2,:));
end

function [Y]=sortidesc(X,countn,ee)
Y = zeros(size(X));
for ii=1:countn-1
    for jj=1:countn-1
        if real(X(jj,1)) < real(X(jj+1,1))
            X=swap(X,jj);
        end
        if abs(real(X(jj,1)) - real(X(jj+1,1)))<ee
            if imag(X(jj+1,1)) > imag(X(jj,1))
                X=swap(X,jj);
            end
        end
    end
end
Y(1:2:end,:) = flipud(X(countn/2+1:end,:));
Y(2:2:end,:) = (X(1:countn/2,:));
end