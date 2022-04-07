
tic
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND POINTWISE GROWTH MODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1)
% findst pointwise growth mode closest to lam_start using an inverse power method for the nonlinear operator pencil
[A_ref,p,n,k,lam_start]=spatial_pencil; % loads example, A(:,:,j), j=1:p+1 are coefficients of lambda^(j-1), p is order, n dimension, k dimension of E^u, lam_start the refefrence point from where to search for closest PGM


 % transform the pencil so lam_start is at the origin
A=poly_shift(A_ref,p,lam_start);

%%%%%%%%%%%%%%%%%%%%%%
% numerical parameters
%%%%%%%%%%%%%%%%%%%%%%
M=20; % highest order of eigenspace to be computed in initial approximation, determines polynomial approximation order of operator pencil
% high M gives high confidence that the closest PGM is determined, runs into issues 
%   with roundoff andunderflow though; also most expensive part... 50-100 are reasonable large numbers
% low M gives usually faster convergence 
M_fine=5; % not so high order iteration when already near singularity; 5-15 is good

tol=1e-8; % tolerance in initial approximation; accuracy of PGM; aim for 1e-6 if not sure, should use M iterations
tol_fine=1e-8; %if confident, 1e-8, otherwise 1e-10
tol_subspace=1e-8; % similar to tol_fine, can cause trouble with gmres but ignore
tol_newton=1e-12; % final tolerance after Newton; 1e-12 should be achievable

tries_max=1; % number of tries with random u0 to achieve accuracy in initial iteration; set to 1 unless real trouble
max_fine=6; % maximum number of iterations in fine loop with order of penciil truncated to M_fine
max_newton = 20; % max iterations in Newton; sert to zero if want to avoid Newton

restart_taylor=false;   % if true, restart iteration with subspace computed by Taylor approximation; Taylor alias "true" is faster, 
                        % otherwise use subspace computed from homotopy eu_hom
verbose=1; 

s_f_init=1; % initial step should not be too large when first exploring, problem when large dlam with converging power series, can be close to 1 later but risk hitting abs spec with ambiguity of eigenspaces; larger M allows for larger s_f_init            
s_f_fine=.85; % step fraction in power method, avoid overshoot so chosen close to but below 1

u0=randn(n,1); % starting vector for inverse power iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% define iota starting with reference eigenspaces and their Taylor expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A0=A(:,:,1); % matrix at start

% here start with the unstable eigenspace, given by a half-dimension so n=2*k
% could also give a starting space that is invariant but a mixture of eigenvalues not separated from the 
% others by real part
[U,T]=schur(full(A0),'complex'); % schur form; 
T_0=maxk(real(diag(T)),k+1); 
med_T=(T_0(k+1)+T_0(k))/2;% find k'th and k+1'st  largest elements of diag(T) and average of the two
[U,T] = ordschur(U,T,real(diag(T))>med_T); % order by grouping eigenvalues with most unstable directions first
Eu=eu_taylor_polynomial_pencil(A,p,U,n,k,M,tol);
%%% also compute the stable eigenspace

[U,T]=schur(-full(A0),'complex'); % schur form; 
[U,T] = ordschur(U,T,real(diag(T))>-med_T); % order by grouping eigenvalues with most unstable directions first
Es=eu_taylor_polynomial_pencil(-A,p,U,n,n-k,M,tol);

iota  = [Eu , Es];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial iteration up to order M, key work done here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flag=0;tries=0;
while (flag==0) && (tries<tries_max)
    [dlam,u_eig,err,M_final,flag,lam_list]=pencil_eig(iota,u0,n,tol,M);
    u0=randn(n,1); % starting vector for inverse power iteration
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
lam_list=lam_list+lam_start;
lam_list_initial=lam_list;
%%%%%%%%%%%%%%%%%%%%%%%%
if true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now do fine iteration with restarts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dlam_o=1e8; % initialize eigenvalue increment for convergence criterion
    it=0;       % initialize counter for number of loops
    while (abs(dlam_o)+abs(dlam)>tol_fine) && (it<max_fine)
        it=it+1;
        dlam_o=dlam;
        if restart_taylor
            % compute approximation for stable/unstable subspaces at new lam location using Taylor expansion, Horner's method
            Meff=(M_fine*(it>1)+M*(it==1));
            dEu=Eu(:,:,Meff+1);
            dEs=Eu(:,:,Meff+1);
            for j=Meff-1:-1:0
                dEu=dEu*dlam+Eu(:,:,j+1);
                dEs=dEs*dlam+Es(:,:,j+1);
            end
            % find ONB  representation and compelements
            [UEu,~]=qr(dEu(:,:,1));
            [UEs,~]=qr(dEs(:,:,1));
            A=poly_shift(A_ref,p,lam); % transform the pencil so lam is at the origin
        
            A0=A(:,:,1);
            UEu=eu_newton(A0,UEu,n,k,tol_subspace);
            UEs=eu_newton(A0,UEs,n,n-k,tol_subspace);
        else
            [UEu,~]=qr(Eu(:,:,1));
            [UEs,~]=qr(Es(:,:,1));
            [UEu,hom_success]=eu_hom(A,p,0,dlam,UEu,n,k,tol_subspace);
            if not(hom_success)
                display(['homotopy of unstable eigenspace did not succeed'])
            end
            [UEs,hom_success]=eu_hom(A,p,0,dlam,UEs,n,n-k,tol_subspace);
            if not(hom_success)
                display(['homotopy of stable eigenspace did not succeed'])
            end
            A=poly_shift(A_ref,p,lam); % transform the pencil so lam is at the origin
        end
        
        % compute the Taylor expansion of this new subspace pair
        Eu=eu_taylor_polynomial_pencil(A,p,UEu,n,k,M_fine,tol_subspace);
        Es=eu_taylor_polynomial_pencil(-A,p,UEs,n,n-k,M_fine,tol_subspace);

        % form iota and it's Taylor jet and do inverse power iteration
        iota  = [Eu , Es];
        [~,Si,~]=svd(iota(:,:,1));
        if min(diag(Si))<1e-12 % already have kernel of iota
    %          break
        end
        [dlam,u_eig,err,M_final,flag,lam_list_n]=pencil_eig(iota,u0,n,tol_fine,M_fine);
        % update current guess for eigenvalue
        size(lam_list_n)
        dlam=s_f_fine*dlam;
        lam_list=[lam_list;lam_list_n+lam];
        lam=lam+dlam;
        if verbose 
            display([num2str(it) '. fine sequence, ' num2str(M_final) ' iters, current lam is ' num2str(lam) ' with last increment ' num2str(dlam)])
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now wrap up with Newton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find pinched double root using Newton, use initial guess for lambda and u_eig to find nu first, then just Newton for nonl evp in nu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    it=0;
    lam0=lam;
    u0=u_eig(:,M_final);

    % compute approximation for stable/unstable subspaces at new lam location using Taylor expansion
    % could use homotopy method here but should be close enough so Taylor works
    dEu=Eu(:,:,M_fine+1);
    dEs=Eu(:,:,M_fine+1);
    for j=M_fine:-1:1
        dEu=dEu*dlam+Eu(:,:,j);
        dEs=dEs*dlam+Es(:,:,j);
    end

    % find ONB  representation and complements
    [UEu,~]=qr(dEu(:,:));
    [UEs,~]=qr(dEs(:,:));

    % confirm approximation with Newton method
    A=poly_shift(A_ref,p,lam0); % transform the pencil so lam is at the origin
    A0=A(:,:,1);
    [Eu,flag]=eu_newton(A0,UEu,n,k,tol_subspace);
    [Es,flag]=eu_newton(A0,UEs,n,n-k,tol_subspace);


    % form iota and it's Taylor jet and do inverse power iteration
    iota0 = [Eu(:,1:k) , Es(:,1:n-k)];

    % now find the intersection of Eu and Es, alias the "almost" kernel of iota
    [Ui,Si,Vi]=svd(iota0);
    k0=find(diag(Si)<1e-2*Si(1,1),1); % here some criterion that identifies the smallest singular values
    if isempty(k0)
        k0=n;
    end
    iotaker=Vi(:,k0:n); % identify the "almost" kernel
    iotaker(1:k,:)=0; % look for the kernal in the stable subspace only, it's replicated only in the unstable subspace
    Vker=iota0*iotaker; % this is the basis of the almost kernel
    [Vker,~,~]=svd(Vker,'econ'); % make Vker orthonormal columns, so projection on kernel is just transpose
    [Ueig,Deig]=eig(Vker'*A0*Vker); % find eigenvectors and eigenvalues of A0 restricted to the kernel of iota, these are the nu candidates

    % now go through the candidates for eigenvectors and nu values and solve the pinched double root equation using Newton
    for sel=n-k0+1:-1:1  % go backwards, starting with the smallest singular value
        lam_newton=[];
        nu0=Deig(sel,sel);
        if verbose 
            display([num2str(lam) '  lambda after fine iteration '])
            display([num2str(nu0) '  Initial guess for nu before Newton iteration ']) 
        end
        e0=Vker*Ueig(:,sel);
    %%%%%%%%%%%%%%%%%
        % now define newton objective function
        %%%%%%%%%%%%%%%%%
        F=@(u,v,lam,nu) [pol_eval(A_ref,p,lam)*u-nu*u;...       % eigenvector to nu of pencil
                            pol_eval(A_ref,p,lam)*v-nu*v-u;...  % Principal vector to eigenvalue
                            e0'*u-1;...                         % normalize eigenvector
                            e0'*v];                             % normalize principal vector
        DF=@(u,v,lam,nu) ...
                [pol_eval(A_ref,p,lam)-nu*eye(n), zeros(n,n) , pol_der_eval(A_ref,p,lam)*u , -u;...
                    -eye(n)  , pol_eval(A_ref,p,lam)-nu*eye(n)  ,  pol_der_eval(A_ref,p,lam)*v , -v;...
                    e0' , zeros(1,n+2);...                     
                    zeros(1,n),e0',0,0];                        
                    % derivative is explicit, so use it; note eveything is complex linear/analytic
        % now defint function formally as function of n+2-vector
        F0=@(z)   F(z(1:n),z(n+1:2*n),z(2*n+1),z(2*n+2));
        DF0=@(z) DF(z(1:n),z(n+1:2*n),z(2*n+1),z(2*n+2));
        % define initial guess as n+2-vector
        z0=[e0;zeros(n,1);lam0;nu0];
        %%% now come the standard Newton iterations %%%%%%%%
        % start of the Newton iterations, compute residual
        res=F0(z0);
        j=0;
        while (norm(res)>tol_newton)&&(j<max_newton)
            j=j+1;
            z0=z0 -DF0(z0)\res;
            res=F0(z0);
            lam_newton=[lam_newton;z0(2*n+1)];
        end
        % collect the relevant information
        lam0=z0(end-1);
        nu0=z0(end);
        if norm(res)>tol_newton
            display(['res is ' num2str(norm(res)) ' above tolerance ' num2str(tol_newton) '!!'])
        else
            display([ num2str(lam0) '  lambda after Newton iteration '])
            display([ num2str(nu0) '  nu after Newton iteration '])
        end
    end
end
lam_list_all=[lam_list;lam_newton];

name=['data_plotting/sh_lam_start_' num2str(lam_start) '_M_' num2str(M_fine) '.mat'];
save(name,'lam_start','lam_list','lam_list_all','lam_list_initial')
plot(log10(abs([lam_list;lam_newton])),'.-')
hold on
plot(log10(abs(lam_list)),'.-')
hold off
toc

