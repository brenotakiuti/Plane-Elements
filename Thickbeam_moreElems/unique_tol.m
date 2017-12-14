function oc = unique_tol(A,tol)
% A MUST be a vector
% Checks how many times a given number appears inside a vector by checking
% if the difference between two positions is smaller than toli.
% A = [1 1.0002 1 2.0002 1i 2i 1.002i 1.1i];
oc = zeros (1,length(A));
ee = 1e-6;
for ii=1:length(A)
    focus = A(ii);
    toli = abs(focus*tol);
    for jj=1:length(A)
        if abs(real(focus))>ee && abs(imag(focus))<ee
            if real(A(jj))*real(focus)>0
                if abs(real(A(jj)) - real(focus)) <toli
                    oc(jj) = oc(jj) + 1;
                end
            end
        elseif abs(imag(focus))>ee && abs(real(focus))<ee
            if imag(A(jj))*imag(focus)>0
                % check how to make s(15)== s(18)
                if abs(imag(A(jj)) - imag(focus)) <toli
                    oc(jj) = oc(jj) + 1;
                end
            end
        else
            if imag(A(jj))*imag(focus)>0
                if abs(imag(A(jj)) - imag(focus)) <toli &&  abs(real(A(jj)) - real(focus)) <toli
                    oc(jj) = oc(jj) + 1;
                end
            end
        end
    end
end
