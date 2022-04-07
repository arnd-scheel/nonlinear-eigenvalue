function [U,flag,count]=eu_hom(A,p,lam1,lam2,U,n,k,tol)

    % continue invariant subspace from lam1 to lam2 for polynomial symbol A; make generic arc perturbation of straight line from lam1 to lam2

    % INPUT
    % matrix polynomial A, size n, dimension of subspace k, tolerance tol, U orthogonal matrix containing
    % approximnate subspace in first k columns and orthocomplement in remaining columns, at point lam1, 

    % OUTPUT
    % E: orthogonal matrix E with invariant subspace in first k columns and orthocomplement in rest
    % at parameter lam2
    tau=0;
    dtau=1e-2;
    dtau_old=1e-2;
    deform=5;(rand-.5); % how much of an arc, random deformation
    % current updated point on arc
    count=0;
    dU=0*U;
    s=@(tau) lam1+tau*(lam2-lam1)+1i*tau.*(1-tau).*(lam2-lam1)*deform;
    while tau+dtau<1
        tau_n=tau+dtau;
        A0=pol_eval(A,p,s(tau_n));
        U_init=U+(s(tau_n)-s(tau))*dU;
        [U_new,newt_success,steps,step_init]=eu_newton(A0,U_init,n,k,tol,false);
        count=count+1;
        if steps>2 || norm(step_init)>1e-3
            dtau=dtau/2;
        else 
            dU=(U_new-U)/(s(tau_n)-s(tau));
            U=U_new;
            tau=tau_n;
        end
    end
    tau=1;
    s=lam2;
    A0=pol_eval(A,p,s);
    [U,newt_success,steps,res_init]=eu_newton(A0,U,n,k,tol,false);
    flag=newt_success;
%      U(:,2)'*A0*U(:,1)
return

