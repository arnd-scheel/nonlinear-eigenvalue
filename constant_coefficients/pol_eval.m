function A_val = pol_eval(A,p,lam)
% evaluates a polynomial, needed in the Newton routine in pdr.m to define objective function inline
% just a simple Horner scheme
A_val=A(:,:,p+1);
    for j=p:-1:1
        A_val=A_val*lam+A(:,:,j);
    end
return
