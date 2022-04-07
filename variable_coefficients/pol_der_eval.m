function A_val=pol_der_eval(A,p,lam)
% evaluates derivative of polynomial, needed in the Newton routine in pdr.m to define derivative of objective function inline
% just a simple Horner scheme
S=size(A);
n=S(1);
if S(2)~=(p+1)*n
    display('errpr in size')
end

A_val=p*A(1:n,(1:n)+p*n);
for j=p:-1:2
    A_val=A_val*lam+(j-1)*A(1:n,(1:n)+(j-1)*n);
end
return
