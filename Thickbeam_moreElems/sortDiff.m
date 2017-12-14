function [y] = sortDiff(y,x,tol)
% y must be a struct with:
%
% PhiQp PhiQn PhiFp PhiFn PsiQp PsiFp PsiQn PsiFn lp ln s kp kn
%
% Sort the matrix k by calculating the angle as function of f. If there is
% an angle where a(n)-a(n-1) > tol, the rows of the matrices are swapped by the
% row with the closest angle to the angle a(n)
% The number of rows in y must be the same number of elements in x.

% This is a fine tuned code, so it works in a very specific and restricted
% condition. Of course it can be fine tunned for each case. I tried to
% stablish a value for "tol" and "tol2" in which the code could work in
% any situation, but it was not tested. The tol is the most important
% condition. It looks at the diff matrix and checks which rows (or
% wavemodes) presents a discontinuity which basically means a bad sorting
% of the k (wavenumbers or eigenvalues) array. The if statement detects a
% discontinuity by comparing one value with its next value in the freq
% dimension of the matrix. If the difference is bigger than tol, it is
% considered a dicontinuity. This difference is in general (when there is
% no discontinuity) defined by the size of the x array in "radSlope". The
% smaller the deltax, the bigger are the differences bewteen one freq and
% the next. By tweeking this x array, I tried to make it so that the
% differences are comparativelly bigger when there is a discontinuity, than
% when there is simply a smooth increase of the slope value (which is
% normal and expercted).

if nargin<3
tol = 1e-1;
end
tol2 = 1e0;
%% Positive wavenumbers
[nr,nc] = size(y.kp);

diff = zeros(nr,nc);


% diff
% simmze(diff)
start = y.PureN+1;
for ii = start:nr
    %     ii
    for jj = start:nr-1
        %         jj
        for mm=start:nr
            diff(mm,:) = radSlope(x,real(y.kp(mm,:)));
        end
        for kk = 1:nc-1
            if abs(diff(jj,kk)-diff(jj,kk+1))>tol && diff(jj,kk)>1e-8
                store = real(y.kp(jj,kk));
                for ll = jj+2:nr-1
                    if abs(real(y.kp(ll,kk+2))-store)<tol2
                        order = 1:1:nr;
                        %                         jj
                        %                         ll
                        [diff(:,kk+1:end),order] = swapMultiple(diff(:,kk+1:end),order,jj,ll);
                        y.kp(:,kk+1:end) = y.kp(order,kk+1:end);
                        y.PhiQp(:,:,kk+1:end) = y.PhiQp(:,order,kk+1:end);
                        y.PhiFp(:,:,kk+1:end) = y.PhiFp(:,order,kk+1:end);
                        y.PsiQp(:,:,kk+1:end) = y.PsiQp(order,:,kk+1:end);
                        y.PsiFp(:,:,kk+1:end) = y.PsiFp(order,:,kk+1:end);
                        y.lp(:,kk+1:end) = y.lp(order,kk+1:end);
                        break;
                    end
                end
                break;
            end
        end
    end
end
%% Negativo wavenumbers
[nr,nc] = size(y.kn);

diff = zeros(nr,nc);

for ii = start:nr
    %     ii
    for jj = start:nr-1
        %         jj
        for mm=start:nr
            diff(mm,:) = radSlope(x,real(y.kn(mm,:)));
        end
        for kk = 1:nc-1
            if abs(diff(jj,kk)-diff(jj,kk+1))>tol && diff(jj,kk)>1e-8
                store = real(y.kn(jj,kk));
                for ll = jj+2:nr-1
                    if abs(real(y.kn(ll,kk+2))-store)<tol2
                        order = 1:1:nr;
                        %                         jj
                        %                         ll
                        [diff(:,kk+1:end),order] = swapMultiple(diff(:,kk+1:end),order,jj,ll);
                        y.kn(:,kk+1:end) = y.kn(order,kk+1:end);
                        y.PhiQn(:,:,kk+1:end) = y.PhiQn(:,order,kk+1:end);
                        y.PhiFn(:,:,kk+1:end) = y.PhiFn(:,order,kk+1:end);
                        y.PsiQn(:,:,kk+1:end) = y.PsiQn(order,:,kk+1:end);
                        y.PsiFn(:,:,kk+1:end) = y.PsiFn(order,:,kk+1:end);
                        y.ln(:,kk+1:end) = y.ln(order,kk+1:end);
                        break;
                    end
                end
                break;
            end
        end
    end
end
end

