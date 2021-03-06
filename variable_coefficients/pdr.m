clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND POINTWISE GROWTH MODES ----- x-dependent coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1) % for reproducability initialize rand
% find pointwise growth mode closest to lam_start using an inverse power method for the nonlinear operator pencil
tic
[A_ref0,Ap_ref,Am_ref,p,n,k,N,x,dx,lam_start]=spatial_pencil; %loads example

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build lambda-dependent equations in bulk of domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% form derivative matrix
e=ones(n*N,1);
D=spdiags([1/24*e -9/8*e 9/8*e -1/24*e],-1:2,N,N+1);
D(1,1:4)=D(1,1:4)+(1/24)*[4 -6 4 -1];
D(N,N-2:N+1)=D(N,N-2:N+1)+(-1/24)*[ -1 4 -6 4]; 
D=D/dx;

F=1/16*spdiags([-e 9*e 9*e -e],-1:2,N,N+1);
F(1,1:4)=F(1,1:4)-(1/16)*[4 -6 4 -1];
F(N,N-2:N+1)=F(N,N-2:N+1)-1/16*[-1 4 -6 4];

D=kron(D,speye(n,n));
F=kron(F,speye(n,n));


% form vector field matrix, average values
A_ref=sparse(n*N,n*(N+1)*(p+1));
A_ref=F*A_ref0;
A_ref(1:n*N,1:n*(N+1))=A_ref(1:n*N,1:n*(N+1))-D;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shift to reference point, A, Ap, and Am
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=poly_shift(A_ref,n*N,n*(N+1),p,lam_start);
Ap=poly_shift(Ap_ref,n,n,p,lam_start);
Am=poly_shift(Am_ref,n,n,p,lam_start);

%%%%%%%%%%%%%%%%%%%%%%
% numerical parameters
%%%%%%%%%%%%%%%%%%%%%%

M=50 ; % highest order of eigenspace to be computed in initial approximation, determines polynomial approximation order of operator pencil
% high M gives high confidence that the closest PGM is determined, runs into issues 
%   with roundoff andunderflow though; also most expensive part...
% low M gives usually faster convergence
M_fine=10; % not so high order iteration when already near singularity but sufficiently high when finding branch point

tol=1e-6; % tolerance in initial approximation; accuracy of PGM; aim for 1e-6 if not sure, should use M iterations
tol_fine=1e-10; %if confident, 1e-8, otherwise 1e-10
tol_subspace=1e-12; % similar to tol_fine, can cause trouble with gmres but ignore
tol_newton=1e-12; % final tolerance after Newton; 1e-12 should be achievable
kron_form=false; % use kronecker matrix form and direct solver, sometimes more accurate but slow vs gmres representation in eu_taylor and eu_newton

tries_max=1; % number of tries with random u0 to achieve accuracy in initial iteration; set to 1 unless real trouble
max_fine=7; % maximum number of iterations in fine loop with order of penciil truncated to M_fine
max_newton = 20; % max iterations in Newton; set to zero if want to avoid Newton


verbose=1; 

s_f_init=.95; % initial step should not be too large when first exploring, problem when large dlam with converging power series, can be close to 1 later but risk hitting abs spec with ambiguity of eigenspaces; larger M allows for larger s_f_init            
s_f_fine=.9; % step fraction in power method, avoid overshoot so chosen close to but below 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build nonlinear (in lambda) boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ap0=Ap(1:n,1:n); % matrix at start

% here start with the unstable eigenspace, given by a half-dimension so n=2*k
% could also give a starting space that is invariant but a mixture of eigenvalues not separated from the 
% others by real part
[U,T]=schur(full(Ap0),'complex'); % schur form; 
T_0=maxk(real(diag(T)),k+1); 
med_T=(T_0(k+1)+T_0(k))/2;% find k'th and k+1'st  largest elements of diag(T) and average of the two

[U,T]=schur(-full(Ap0),'complex'); % schur form; 
[U,T] = ordschur(U,T,real(diag(T))>-med_T); % order by grouping eigenvalues with most unstable directions first
Es=eu_taylor_polynomial_pencil(-Ap,p,U,n,n-k,M,tol);

Am0=Am(1:n,1:n); % matrix at start

% here start with the unstable eigenspace, given by a half-dimension so n=2*k
% could also give a starting space that is invariant but a mixture of eigenvalues not separated from the 
% others by real part
[U,T]=schur(full(Am0),'complex'); % schur form; 
T_0=maxk(real(diag(T)),k+1); 
med_T=(T_0(k+1)+T_0(k))/2;% find k'th and k+1'st  largest elements of diag(T) and average of the two


[U,T] = ordschur(U,T,real(diag(T))>med_T); % order by grouping eigenvalues with most unstable directions first
Eu=eu_taylor_polynomial_pencil(Am,p,U,n,k,M,tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% put pencil together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% makes a square matrix pencil of size n*(N+2) and order M (that is entries up to M+1
for j=1:p+1
    iota(1:n*(N+2),(1:n*(N+2))+(j-1)*n*(N+2))=   [A(1:n*N,(1:n*(N+1))+(j-1)*n*(N+1)), sparse(n*N,n);...
                    (j==1)*speye(n,n), sparse(n,n*N), sparse(n,n-k),sparse(Eu(1:n,(1:k)+(j-1)*(k)));...
                    sparse(n,n*N),(j==1)*speye(n,n), sparse(Es(1:n,(1:n-k)+(j-1)*(n-k))), sparse(n,k)];
end
for j=p+2:M+1
     iota(1:n*(N+2),(1:n*(N+2))+(j-1)*n*(N+2))=   [sparse(n*N,n*(N+2));...
                    speye(n,n), sparse(n,n*N), sparse(n,n-k),sparse(Eu(1:n,(1:k)+(j-1)*(k)));...
                    sparse(n,n*N),speye(n,n), sparse(Es(1:n,(1:n-k)+(j-1)*(n-k))), sparse(n,k)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% define iota starting with reference eigenspaces and their Taylor expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u0=randn(n*(N+2),1); % starting vector for inverse power iteration

u0=ones(n*(N+2),1);
% here start with the unstable eigenspace, given by a half-dimension so n=2*k
% could also give a starting space that is invariant but a mixture of eigenvalues not separated from the 
% others by real part

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial iteration up to order M, key work done here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
flag=0;tries=0;
while (flag==0) && (tries<tries_max)
    [dlam,u_eig,err,M_final,flag,lam_list]=pencil_eig(iota,u0,n*(N+2),tol,M);
    lam_list=lam_list+lam_start;
    u0=randn(n*(N+2),1);  % starting vector for inverse power iteration
    tries=tries+1;
end
if flag==0 && verbose
    display(['error too large = ' num2str(err) ' after initial iteration, try larger M or smaller tol'])
%      return
end
dlam=s_f_init*dlam; % slightly reduce step size
lam=lam_start+dlam;% this one does better with a slightly smaller step to avoid overshooting and getting problems with non-con ergent E^j
if verbose 
    display(['initially ' num2str(M_final) ' iterations out of ' num2str(M)])
    display([ num2str(lam)  '  lambda after initial iteration ']) 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now do fine iteration with restart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dlam_o=1e8; % initialize eigenvalue increment for convergence criterion
it=0; % initialize counter for number of loops
while (abs(dlam_o)+abs(dlam)>tol_fine) && (it<max_fine)
    u0=u_eig(:,1);
%      u0=randn(n*(N+2),1);
%      u0=ones(n*(N+2),1);
    it=it+1;
    dlam_o=dlam;
    
    A=poly_shift(A_ref,n*N,n*(N+1),p,lam);
    Ap=poly_shift(Ap_ref,n,n,p,lam);
    Am=poly_shift(Am_ref,n,n,p,lam);
    
    % compute approximation for stable/unstable subspaces at new lam location using Taylor expansion, Horner's method
    Meff=(M_fine*(it>1)+M*(it==1));
    dEu=Eu(1:n,(1:k)+k*Meff);
    dEs=Es(1:n,(1:(n-k))+(n-k)*Meff);
    for j=Meff-1:-1:0
        dEu=dEu*dlam+Eu(1:n,(1:k)+k*j);
        dEs=dEs*dlam+Es(1:n,(1:(n-k))+(n-k)*j);
    end
    
    % find ONB  representation and complements
    [UEu,~]=qr(Eu(1:n,1:k));
    [UEs,~]=qr(Es(1:n,1:n-k));
%      A0=A(:,:,1);
    [UEu,flagu]=eu_newton(Am(1:n,1:n),UEu,n,k,tol_subspace);
    [UEs,flags]=eu_newton(Ap(1:n,1:n),UEs,n,n-k,tol_subspace);
    if not(flagu && flags)
        display('restart failed since subspaces did not converge: try doing subspace homotopy or smaller step s_init/s_fine instead')
        break
    end
    % compute the Taylor expansion of this new subspace pair
    Eu=eu_taylor_polynomial_pencil(Am,p,UEu,n,k,M_fine,tol_subspace);
    Es=eu_taylor_polynomial_pencil(-Ap,p,UEs,n,n-k,M_fine,tol_subspace);
    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % put pencil together
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % makes a square matrix pencil of size n*(N+2) and order M (that is entries up to M+1
    for j=1:p+1
        iota(1:n*(N+2),(1:n*(N+2))+(j-1)*n*(N+2))=   [A(1:n*N,(1:n*(N+1))+(j-1)*n*(N+1)), sparse(n*N,n);...
                        (j==1)*speye(n,n), sparse(n,n*N), sparse(n,n-k),sparse(Eu(1:n,(1:k)+(j-1)*(k)));...
                        sparse(n,n*N),(j==1)*speye(n,n), sparse(Es(1:n,(1:n-k)+(j-1)*(n-k))), sparse(n,k)];
    end
    for j=p+2:M_fine+1
        iota(1:n*(N+2),(1:n*(N+2))+(j-1)*n*(N+2))=   [sparse(n*N,n*(N+2));...
                        speye(n,n), sparse(n,n*N), sparse(n,n-k),sparse(Eu(1:n,(1:k)+(j-1)*(k)));...
                        sparse(n,n*N),speye(n,n), sparse(Es(1:n,(1:n-k)+(j-1)*(n-k))), sparse(n,k)];
    end

% makes a square matrix pencil of size n*(N+2) and order M (that is entries up to M+1
%      [~,Si,~]=svd(full(iota(1:n*(N+2),1:n*(N+2))));
%      if min(diag(Si))<1e-12 % already have kernel of iota
%          break
%      end
    [dlam,u_eig,err,M_final,flag,lam_list_n]=pencil_eig(iota,u0,n*(N+2),tol_fine,M_fine);
    lam_list=[lam_list;lam_list_n+lam];
    % update current guess for eigenvalue
    dlam=s_f_fine*dlam;
    lam=lam+dlam;
    if verbose
        display(['Current increment ' num2str(dlam) '  to current value ' num2str(lam)])
    end
end
if verbose 
    display([num2str(it) '  fine iterations out of max ' num2str(max_fine) ' with last increment ' num2str(dlam)])
end

%  plot(lam_list+lam_start,'.-')
plot(real(log10(lam_list+2)),'.-')
%  plot(real((lam_list)),'.-')

%  plot(log10(1:length(lam_list))',real(log10(lam_list+lam_start+2)),'.-')

name=['ac_lam_start_restart_n' num2str(-lam_start) '_M_' num2str(M) '.mat'];

save(name,'lam_start','lam_list','M','M_fine')
toc

return




toc
return




