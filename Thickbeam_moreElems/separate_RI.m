function [R_part, I_part] = separate_RI(X,tol)

ee = tol;
n = length(X);

count_R = 1;
count_I = 1;

R_part = 0;
I_part = 0;

for ii=1:n
    if abs(real(X(1)))>ee && abs(imag(X(1)))<ee
        % Is REAL
        R_part(count_R) = X(1);
        count_R = count_R+1;
        X(1) = [];
    elseif abs(imag(X(1)))>ee && abs(real(X(1)))<ee
        % Is IMAGINARY
        I_part(count_I) = X(1);
        count_I = count_I+1;
        X(1) = [];
    else
        error(['Complex value found, in position [' num2str(ii) '] either eliminate the complex values or loosen the tolerance']);
    end
end

end