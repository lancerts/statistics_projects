function [num, supp] = Func_CalcNz(X, offDiag, tol)

if offDiag == 0
    val = X;
else
    val = X - diag(diag(X));
end
supp = find(abs(val)>= tol);
num = length(supp);