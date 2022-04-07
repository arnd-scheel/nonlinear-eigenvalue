

function [dlam,u_eig,err,M_final,flag,lam_list]=pencil_eig(iota,u0,n,tol,M)
    % simlar to pencil_eig_2 or analytic root finder
    % INPUT
    % iota: pencil, iota(lam)=iota(:,:,1) + iota(:,:,2)*lam+ iota(:,:,3)*lam^2 +... + iota(:,:,M)*lam^(M-1) matrix at order lambda=0
    % u0: starting vector nx1
    % n: size of problem n, 
    % tol: tolerance, 
    % M: max number of iterations, that is, order to which iota_n is provided
    % OUTPUT
    % dlam: final guess for current eigenvalue
    % u_eig:    associated eigenvector
    % err: last increment in eigenvector and lack of colinearity in power method
    % M_final: at most M, but less if converged to desired tolerance earlier
    % flag: true if converged

    u=zeros(n,M); % vector of iterates
    u(:,1)=u0;       % initialize first value
    lam_list=[];

    gam_o=1;    % old estimate for root
    gam_n=0;     % new estimate for root

    %%%%%%%%%%%     MAIN ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j=2:M
        gam_o=gam_n;
        r=zeros(n,1);
        for jj=j-1:-1:1 
            r=r+iota(:,:,jj+1)*u(:,jj);
        end
        u(:,1:j) = [-iota(:,:,1)\r, u(:,1:j-1)];% define new iterate
        u=u/norm(u,Inf);  % renormalize to avoid blowup of u
        if isnan(u)
            break
        end
        gam_n=(u(:,1)'*u(:,2))/(u(:,1)'*u(:,1));% find eigenvalue candidate
        lam_list=[lam_list;gam_n];
        err_gam=abs(gam_n-gam_o);
        err_eig=norm(gam_n*u(:,1)-u(:,2));%/norm(-u(:,j-1));
%          if (j>3) && err_gam+err_eig<tol % at least 4 iterations
%              break
%          end
    end 
    dlam=gam_n;
    u_eig=u(:,1:j)/norm(u(:,j),'inf');
    M_final=j;
    
    err=err_eig+err_gam;
    flag=err<tol;
return
    
