function [M] = cleanComplex(M,ee)

[r,c] = size(M);

for jj = 1:c
    for ii = 1:r
        if abs(real(M(ii,jj)))<0+ee
            M(ii,jj) = imag(M(ii,jj))*1i;
        end
        if abs(imag(M(ii,jj)))<0+ee
            M(ii,jj) = real(M(ii,jj));
        end
    end
end