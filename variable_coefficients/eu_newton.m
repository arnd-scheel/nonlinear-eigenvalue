

function [E,flag,steps,step_init]=eu_newton(A,U,n,k,tol,kron_form)

% find invariant subspace near approximate invariant subspace
% INPUT
% matrix A, size n, dimension of subspace k, tolerance tol, U orthogonal matrix containing
% approximnate subspace in first k columns and orthocomplement in remaining columns

% OUTPUT
% E: orthogonal matrix E with invariant subspace in first k columns and orthocomplement in rest
% flag: true if convergence after 20 iterations, false otherwise
    At=U'*A*U; % change coordinates so subspace is canonical hyperspace of first k coordinates

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now compute expansion

    % define matrix blocks and homological equation for invariance, with subspace given as graph over first k coordinates, so solve for matrix H R^k-> R^n-k

    A00 = At(1:k,1:k);
    A01 = At(1:k,k+1:n);
    A10 = At(k+1:n,1:k);
    A11 = At(k+1:n,k+1:n);
    % invariance condition for graph H as g=0
    f=@(H) A10+A11*H-H*A00-H*A01*H;
    
    % now initialize Newton
    H0=zeros(n-k,k);
    res=f(H0);
    j=0;
    step_init=0;
    % standard Newton routine
    while (norm(res,'inf')>tol) && (j<20)
        % linearization, lots of reshaping to represent matrix multiplication in column vectors
        dH=sylvester(A11-H0*A01,-(A00+(A01*H0)),res);
        if j==0
            step_init=norm(dH);
        end
        H0=H0-dH; % implement Newton step
        res=f(H0); % compute new residual
        j=j+1;
    end
    flag = norm(res,'inf')<tol; 
    % write actual subspace as [id;H0], then change back coordinates, and find ONB E
    [E,~]=qr(U*[eye(k,k);H0]);
    steps=j;
return



return
