function A_val = pol_eval(A,p,lam)
% evaluates a polynomial, needed in the Newton routine in pdr.m to define objective function inline
% just a simple Horner scheme
% sparse version
S=size(A);
n=S(1);
if S(2)~=(p+1)*n
    display('errpr in size')
end

A_val=A(1:n,(1:n)+p*n);
for j=p:-1:1
    A_val=A_val*lam+A(1:n,(1:n)+(j-1)*n);
end
return
