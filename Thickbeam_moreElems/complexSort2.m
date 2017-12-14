function [Pos,Neg]=complexSort2(X,MODE,tol)
% Works symilarly to the command sort but sorts complex numbers considering
% the real part as always bigger than the imaginary. Ex.:
%               1+1i>0+2i
%               1+1i>0+1.4142i
%               1+0i>0+1000i
% If the 'descent'  mode is chosen, the numbers with bigger real parts will
% always be in the first half of the sorted vector.
% Only works for column vectors
%
% this version also checks if the absolute value also follows the sorting,
% the vectors must be ordered both as a complex AND the absolute value.
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
ee = 1e-3;

if strcmp(MODE,'descend') || strcmp(MODE,'Descend')
    for ii=1:n-1
        for jj=1:n-1
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
    
    Pos = round(X(1:n/2,:));
    Negi = round(X(n/2+1:end,:));
%     size(Pos)
%     size(Negi)
    
    
    [nNeg,~] = size(Negi);
    if nNeg==1
        Neg = Negi;
    elseif rem(nNeg,2) ~= 0
            Pos(end,:)=[];
            Negi(1,:)=[];
        else
        for ii = 0:2:nNeg-1
            ii;
            Neg(ii+1:ii+2,:) = Negi(nNeg-ii-1:nNeg-ii,:);
        end
    end 
    
%     for jj=1:nNeg-1
% %                 jj
% %                 Neg
% %                 Pos
% %                 abs(Neg)
% %                 abs(Neg(jj,1))
% %                 abs(Neg(jj+1,1))
% %                 abs(Neg(jj,1)) > abs(Neg(jj+1,1))
% %                 round(abs(Neg(jj,1))) > round(abs(Neg(jj+1,1)))
%         if round(abs(Neg(jj,1))) < round(abs(Neg(jj+1,1)))
%             Neg=swap2(Neg,jj);
%             Pos=swap2(Pos,jj);
%         end
%     end
%     Neg
%     Pos
    
elseif strcmp(MODE,'ascend') || strcmp(MODE,'Ascend')
    
    for ii=1:n-1
        for jj=1:n-1
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
    
    Pos = (X(n/2+1:end,:));
    Negi = (X(1:n/2,:));
%     size(Pos)
%     size(Negi)

    [nNeg,~] = size(Negi);
    if rem(nNeg,2) ~= 0
        Pos(1,:)=[];
        Negi(end,:)=[];
    end
    
    [nNeg,~] = size(Negi);
    if nNeg==1
        Neg = Negi;
    else
        for ii = 0:2:nNeg-1
            ii;
            Neg(ii+1:ii+2,:) = Negi(nNeg-ii-1:nNeg-ii,:);
        end
    end
    
%     for jj=1:nNeg-1
%         %         jj;
%         %         Neg;
%         %         abs(Neg)
%         %         abs(Neg(jj,1))
%         %         abs(Neg(jj+1,1))
%         %         abs(Neg(jj,1)) > abs(Neg(jj+1,1))
%         %         round(abs(Neg(jj,1))) > round(abs(Neg(jj+1,1)))
%         if round(abs(Neg(jj,1))) > round(abs(Neg(jj+1,1)))
%             Neg=swap2(Neg,jj);
%             Pos=swap2(Pos,jj);
%         end
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

function [X]=swap2(X,ii)
swp1 = X(ii-1,:);
swp2 = X(ii,:);
X(ii-1,:)=X(ii+1,:);
X(ii,:)=X(ii+2,:);
X(ii+1,:)=swp1;
X(ii+2,:)=swp2;
% swp = I(ii);
% I(ii)=I(ii+1);
% I(ii+1)=swp;
end