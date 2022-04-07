function [U,flag]=eu_hom(A,p,lam1,lam2,U,n,k,tol)

    % continue invariant subspace from lam1 to lam2 for polynomial symbol A; make generic arc perturbation of straight line from lam1 to lam2

    % INPUT
    % matrix polynomial A, size n, dimension of subspace k, tolerance tol, U orthogonal matrix containing
    % approximnate subspace in first k columns and orthocomplement in remaining columns, at point lam1, 

    % OUTPUT
    % E: orthogonal matrix E with invariant subspace in first k columns and orthocomplement in rest
    % at parameter lam2
    S=25;
    newt_success=false;
    while not(newt_success) && S<=5*2^10
        tau=linspace(0,1,S);
        s=lam1+tau*(lam2-lam1)+1i*tau.*(1-tau).*(lam2-lam1)*(rand-.5);
        j=1;
        newt_success=true;
        while (j<=length(s)) && newt_success 
            A0=pol_eval(A,p,s(j));
            [U,newt_success]=eu_newton(A0,U,n,k,tol,false);
            j=j+1;
        end
        if not(newt_success)
            display(['repeating subspace homotopy with ' num2str(S) ' steps'])
            S=2*S;
        end
    end
    flag=newt_success;
return

