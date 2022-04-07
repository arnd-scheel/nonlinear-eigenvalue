function Eu=eu_taylor_polynomial_pencil(A,p,U,n,k,M,tol)


    %%%%%%%% compute expansion of unstable eigenspace to generalized 
    %   eigenvalue problem depending on parameter lambda 
    %   !!!!!!!!!!!!!!!!!only linear in lambda!!!!!!!!!!!!!!!!!!    
    % see 'Grassmannian_expansion_notes.pdf' for theory
    % 
    % arguments are
    %   A(:,:,j) coefficients of polynomial pencil sum A(:,:,j+1) lam^j  
    %   p order of polynomial, so A is of size n x n x (p+1)
    %   U: orthogonal matrix with invariant subspace in first k columns, its orthocomplement
    %      in the remaining columns; if only approximate subspace known first correct with eu_newton
    %   n: size of square matrices A,B,U
    %   k: dimension of invariant subspace
    %   M: order to which subspace is to be expanded
    %
    %   returns:
    %   Eu0: invariant subspace, simply U(:,1:k)
    %   Eu:  such that Eu(lam)=Eu0+lam*Eu(:,:,1)+lam^2*Eu(:,:,2)+lam^3*Eu(:,:,3)+...

    %  Example input
    %      lams=1i; % near which point do we want expansion
    %      B=zeros(4,4);B(4,1)=1; % how does the (linear) parameter enter
    %      A=[0 1 0 0;0 0 1 0; 0 0 0 1;0 0 0 0]+lams*B; % matrix at lambda=lambda=lambda*
    %  
    %      n=4; % size of matrices
    %      k=2; % dimension of unstable eigenspace
    %       
    %       to obtain U for the unstable eigenspace then do
    %       [U,T]=schur(A,'complex'); % schur form; 
    %       [U,~] = ordschur(U,T,real(diag(T))>median(real(diag(T)))); % order by grouping      
    %       eigenvalues with most unstable directions first

    %  
    %      %%%%
    %      M=6; % highest order of eigenspace to be computed

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    % normalize the problem at lambda=lams, with tilde matrices given in new coordinates defined through 
    % coordinate transformation U
    % define matrix blocks

    for j=0:p
        At(:,:,j+1)=U'*A(:,:,j+1)*U;
        A00(:,:,j+1) = At(1:k,1:k,j+1);
        A01(:,:,j+1) = At(1:k,k+1:n,j+1);
        A10(:,:,j+1) = At(k+1:n,1:k,j+1);
        A11(:,:,j+1) = At(k+1:n,k+1:n,j+1);
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now compute expansion; first order terms are special
    H(:,:,1)=sylvester(A11(:,:,1),-A00(:,:,1),-A10(:,:,2));
    % compute higher order expansions

    for ell=2:M
        % define residual on right-hand side
        if ell<=p 
            R=-A10(:,:,ell+1);
        else 
            R=0*A10(:,:,p+1);
        end
        for m=1:min(p,ell-1)
                R=R+H(:,:,ell-m)*A00(:,:,m+1)-A11(:,:,m+1)*H(:,:,ell-m);
        end
        for m=0:min(p,ell-1)
                for i=1:ell-1-m
                    R=R+H(:,:,i)*A01(:,:,m+1)*H(:,:,ell-i-m);
                end       
        end
        % solve homological equation    
        H(:,:,ell)=sylvester(A11(:,:,1),-A00(:,:,1),R);
    end
    %  transform back to original coordinates
    Eu(:,:,1)=U(:,1:k);
    for ell=1:M
        Eu(:,:,ell+1)=U(:,k+1:n)*H(:,:,ell);
    end
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      % now compute expansion; first order terms are special
%           
%      % form kronecker product
%      [U11,R11]=schur(A11(:,:,1),'complex'); % first do schur to solve sylvester equations 
%      [U00,R00]=schur(A00(:,:,1),'complex'); % in fact schur is better than hess
%      
%      % action of multiplication on column vector through Kronecker product , https://en.wikipedia.org/wiki/Kronecker_product#Kronecker_sum_and_exponentiation bullet matrix equations
%      L_A = kron(speye(k,k),R11)-kron(R00.',speye(n-k,n-k));
%  %          [LL,UU,PP]=lu(L_A);
%      H_vec = L_A\reshape(-U11'*A10(:,:,2)*U00,(n-k)*k,1);
%  %          H_vec = UU\(LL\(PP*reshape(-U11'*B10*U00,(n-k)*k,1)));
%      H=zeros(n-k,k,M);
%      H(:,:,1)=U11*reshape(H_vec,n-k,k)*U00';
%      
%      % compute higher order expansions
%  
%      for ell=2:M
%          % define residual on right-hand side
%          if ell<=p 
%              R=-A10(:,:,ell+1);
%          else 
%              R=0*A10(:,:,p+1);
%          end
%          for m=1:min(p,ell-1)
%                  R=R+H(:,:,ell-m)*A00(:,:,m+1)-A11(:,:,m+1)*H(:,:,ell-m);
%          end
%          for m=0:min(p,ell-1)
%                  for i=1:ell-1-m
%                      R=R+H(:,:,i)*A01(:,:,m+1)*H(:,:,ell-i-m);
%                  end       
%          end
%          % solve homological equation    
%          H_vec = L_A\reshape(U11'*R*U00,(n-k)*k,1);
%          H(:,:,ell)=U11*reshape(H_vec,n-k,k)*U00';
%      
%      end
%  
%      %  transform back to original coordinates
%      Eu(:,:,1)=U(:,1:k);
%      for ell=1:M
%          Eu(:,:,ell+1)=U(:,k+1:n)*H(:,:,ell);
%      end
return


