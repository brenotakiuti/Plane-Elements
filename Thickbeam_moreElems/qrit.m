function [Q,R,count] = qrit(A,tol)

er = 100000;
[Q,R] = qr(A);
count = 0;

while er>tol
    A = R*Q;
    [Qn,Rn] = qr(A);
    er = max(abs(diag(R)-diag(Rn)));
    Q = Qn;
    R = Rn;
    count = count+1;
end

