function A_val=pol_der_eval(A,p,lam)
% evaluates derivative of polynomial, needed in the Newton routine in pdr.m to define derivative of objective function inline
% just a simple Horner scheme
A_val=p*A(:,:,p+1);
    for j=p:-1:2
        A_val=A_val*lam+(j-1)*A(:,:,j);
    end
return
