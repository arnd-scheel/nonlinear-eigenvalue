function [A,Ap,Am,p,n,k,N,x,dx,lam_start]=spatial_pencil;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define examples loaded in pdr.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% equation is u_x=A(x;lam)u, A(x;lam)->Ap/Am(lam) 
% A(x; *), Ap(*), and Am(*) are polynomial of defree p in lam, with expansion as in 

% A(:,:,1) u+lam A(:,:,2) u +...+ lam^p A(:,:,p+1) u=0 
% Ap/Am have k-dimensional unstable subspace at lam_start
% change dx and size of domain here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % %%%  Allen-Cahn linearization, linear in lambda pancil
%  lam_start=-1.8; %with pars M=50, M_fine=10, s_f_init=.95;s_f_fine=.9;
lam_start=-1.755; 
lam_start=5;

%  for resonances choose larger domain
L=7;

N=1400;
%  
% sufficient for normal domain size computations with error 1e-13(!) at 4th order
%  L=8;
%  N=12000;

x=linspace(-L,L,N+1); % N+1 grid opints for x grid, so u will live on N+1 points, N equations
dx=x(2)-x(1);
a=1-3*tanh(x/sqrt(2)).^2; 


%%% if want to eliminate root at branch point  a=a-10*sech(x/.2);

p=1;
n=2;
k=1;

% store tensor Ap(1:n,n,1:p+1) in sparse matrix Ap(:,:,k) goes to Ap(:,:+(k-1)*n)
Ap=sparse(n,n*(p+1));
Ap(1:n,1:n)=sparse([0 1; 2 0]);
Ap(2,1+(2-1)*2)=1; % how does the (linear) parameter enter

% swap signs if wanting to compute resonance   
% Ap=-Ap;

% same asymptotic linearization at pm infinity
Am=Ap;

% Kronecker formulation
A=sparse(n*(N+1),n*(N+1)*(p+1)); 
% store tensor A(1:n*N,n*N,1:p+1) in sparse matrix A(:,:,k) goes to A(:,:+(k-1)*n*N)
for j=1:N+1
    A(2*j-1:2*j,(2*j-1:2*j) +(1-1)*n*(N+1))=sparse([0 1; -a(j) 0]);
    A(2*j-1:2*j,(2*j-1:2*j) +(2-1)*n*(N+1))=sparse([0 0; 1 0]); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % %%%  Allen-Cahn linearization, variations on linear in lambda pencil

% anchor Riemann surface at rho, lam=rho+gamma^2
% ok for rho=0, rho=-1, problems at rho=-2??
%  
%  rho=-2;
%  lam_start_actual=-1.6;
%  lam_start_actual=-1.9;
%  %  lam_start_actual=-.9;
%  
%  lam_start=sqrt(lam_start_actual-rho); %this is really gam_start, recover by computing lam^2-rho
%  
%  L=7;
%  N=500;
%  x=linspace(-L,L,N); % N grid opints for x grid, so u will live on N+1 points, N equations
%  dx=x(2)-x(1);
%  a=1-3*tanh(x/sqrt(2)).^2; 
%  
%  p=2; % 2 for quadratic pencil
%  n=2;
%  k=1;
%  
%  Ap(:,:,1)=[0 1; rho-(-2) 0];
%  Ap(:,:,2)=[0 0; 0 0];
%  Ap(:,:,3)=[0 0; 1 0]; % how does the quadratic parameter enter
%     
%  Am=Ap;
%  
%  % tensor formulation
%  %  for j=1:N
%  %      A(:,:,j,1)=[0 1; -a(j) 0];
%  %      A(:,:,j,2)=zeros(2,2);A(2,1,j,2)=1;
%  %  end
%  
%  % Kronecker formulation
%  A=zeros(n*N,n*N,p+1); % would really want this to be sparse
%  for j=1:N
%      A(2*j-1:2*j,2*j-1:2*j,1)=[0 1; rho-a(j) 0];
%      A(2*j-1:2*j,2*j-1:2*j,2)=[0 0; 0 0];
%      A(2*j-1:2*j,2*j-1:2*j,3)=[0 0; 1 0];
%  end
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%  Allen-Cahn linearization, resolve the branch point
%  lam_start=1.4;
%  L=15;
%  N=150;
%  x=linspace(-L,L,N); % N grid opints for x grid, so u will live on N+1 points, N equations
%  dx=x(2)-x(1);
%  a=1-3*tanh(x/sqrt(2)).^2; 
%  
%  
%  p=2;
%  n=2;
%  k=1;
%  
%  Ap(:,:,1)=[0 1; 0 0];
%  Ap(:,:,2)=zeros(2,2);
%  Ap(:,:,3)=zeros(2,2);
%  Ap(2,1,3)=1; % how does the quadratic parameter enter
%     
%  Am=Ap;
%  
%  
%  A=zeros(n*N,n*N,p+1); % would really want this to be sparse
%  for j=1:N
%      A(2*j-1:2*j,2*j-1:2*j,1)=[0 1; -a(j)+2 0];
%      A(2*j-1:2*j,2*j-1:2*j,2)=[0 0; 0 0];
%      A(2*j-1:2*j,2*j-1:2*j,3)=[0 0; 1 0];
%  end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % %%%  Schroedinger with potential
%  lam_start=0.1+.1*1i;
%  L=10;
%  N=300;
%  x=linspace(-L,L,N); % N grid opints for x grid, so u will live on N+1 points, N equations
%  dx=x(2)-x(1);
%  a=exp(-x.^2); 
%  
%  
%  p=2; % 2 for quadratic pencil
%  n=2;
%  k=1;
%  
%  % anchor Riemann surface at rho, lam=rho+gamma^2
%  % ok for rho=0, rho=-1, problems at rho=-2...??
%  
%  rho=0;
%  
%  Ap(:,:,1)=[0 1; rho 0];
%  Ap(:,:,2)=zeros(2,2);
%  Ap(:,:,3)=zeros(2,2);
%  Ap(2,1,3)=1; % how does the quadratic parameter enter
%     
%  Am=Ap;
%  
%  % tensor formulation
%  %  for j=1:N
%  %      A(:,:,j,1)=[0 1; -a(j) 0];
%  %      A(:,:,j,2)=zeros(2,2);A(2,1,j,2)=1;
%  %  end
%  
%  % Kronecker formulation
%  A=zeros(n*N,n*N,p+1); % would really want this to be sparse
%  for j=1:N
%      A(2*j-1:2*j,2*j-1:2*j,1)=[0 1; rho-a(j) 0];
%      A(2*j-1:2*j,2*j-1:2*j,3)=[0 0; 1 0];
%  end
%  %  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of example definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
