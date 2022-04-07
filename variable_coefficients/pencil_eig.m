

function [dlam,u_eig,err,M_final,flag,lam_list]=pencil_eig(iota,u0,n,tol,M)
    % simlar to pencil_eig_2 or analytic root finder
    % INPUT
    % iota_0: matrix at order lambda=0
    % iota_n: its Taylor series coefficients iota=iota_0 + sum_j=1^M iota_n(:,:,j) lambda^j + ...
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
    iota_0=iota(1:n,1:n);  % to be inverted for iterations in inverse power method
    lam_list=[];
    u=zeros(n,M); % vector of iterates
    u(:,1)=u0;       % initialize first value
    gam_o=1;    % old estimate for root
    gam_n=0;     % new estimate for root
    
    %%%%%%%%%%%     MAIN ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for j=2:M
        gam_o=gam_n;
        r=iota(1:n,n+1:n*j)*reshape(u(1:n,1:j-1),n*(j-1),1);
        u(:,1:j) = [-iota_0\r, u(:,1:j-1)];% define new iterate; most of the time spent here for large n
        u=u/norm(u(:,1),'inf');  % renormalize to avoid blowup of u
        gam_n=(u(:,2)'*u(:,2))/(u(:,2)'*u(:,1));% find eigenvalue candidate
        lam_list=[lam_list;gam_n];
        err_gam=abs(gam_n-gam_o);
        err_eig=norm(gam_n*u(:,1)-u(:,2));%/norm(-u(:,j-1));
        if (j>3) && err_gam+err_eig<tol % at least 4 iterations
            break
        end
    end 
    dlam=gam_n;
    u_eig=u(:,1:j)/norm(u(:,j),'inf');
    M_final=j;
    err=err_eig+err_gam;
    flag=err<tol;
return
    
