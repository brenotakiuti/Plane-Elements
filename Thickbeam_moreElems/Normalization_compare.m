% Run a code and use this to easily normalize the PHI with respect to any
% row of itself, or with respect to the other direction.

[R,C] = size(PhiQa_p1);
for ii = 1:R
    Ap(ii,:) = PhiQa_p1(ii,:)./PhiQa_p1(1,:);
end

for ii = 1:R
    An(ii,:) = PhiQa_n1(ii,:)./PhiQa_n1(1,:);
end

for ii = 1:R
    Bp(ii,:) = PhiFa_p1(ii,:)./PhiFa_p1(1,:);
end

for ii = 1:R
    Bn(ii,:) = PhiFa_n1(ii,:)./PhiFa_n1(1,:);
end