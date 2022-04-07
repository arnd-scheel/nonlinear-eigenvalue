% Finds the first PGM of the Laplacian with Dirichlet boundary conditions; correctly identifies the origin as a PGM although this is not a root of the Evans function in this case

M=300;              % highest derivative prepared, equivalent to highest delay incorporated in the iteration
lam0=.1;            % initial guess for PGM 
n=2;                % dimension of the first-order system
tol=1e-5;           % tolerance
u0=[1;1];           % initial guess for starting vector
 
f=@(lam) sqrt(lam);     % (negative) stable subspace, second component, and Taylor series coefficient derivatives
fn=@(lam,n) 1./factorial(n).*lam.^(1/2-n).*sqrt(pi)/2./gamma(3/2-n); 
                    % the 1/n! is so that f = f0 + f1 x + f2 x^2+ f3 x^3...
                    
iota_0=[0 -1;1 f(lam0)]; 
                    % first column is the Dirichlet subspace, second column the stable eigenspace
                    
% now define the tensor iota needed in the inverse power iteration
for j=1:M
    iota_n(:,:,j)=[0 0;0 fn(lam0,j)];
end

% here is where all the work happens: pencil_eig takes iota0 and its Taylor span, initial guess, dimension, tolerance, and max delay as parameters, returns:
    % dlam, the eigenvalue (relative to anchor point lam0
    % eigenvector resulting from the iteration
    % error (anticipated
    % M_final: how many iterations actually happened? Clearly M_final <= M
    % flag: if true we have that errors are smaller than tolerance (see code for actual definition of errors)
    
[dlam,u_eig,err,M_final,flag]=pencil_eig(iota_0,iota_n,u0,n,tol,M);

M_final             % number of iterations
dlam+lam0           % should be zero!



