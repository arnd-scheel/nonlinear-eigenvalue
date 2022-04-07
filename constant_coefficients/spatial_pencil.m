function [A,p,n,k,lam_start]=spatial_pencil;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define examples loaded in pdr.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% equation is A(:,:,1) u+lam A(:,:,2) u +...+ lam^p A(:,:,p+1) u=0 with k-dimensional unstable subspace at lam_start
% some old example below

% %%% 


%%%% Drift-Diffusion uxx+c ux+ mu u = lam u
%      lam_start=5; % near which point do we want expansion
%      p=1;
%      mu=1;c=2; % double root at mu-c^2/4
%      
%      n=2; % size of matrices
%      k=1; % dimension of unstable eigenspace
%      A(:,:,1)=[0 1; -mu -c];
%      A(:,:,2)=zeros(2,2);A(2,1,2)=1; % how does the (linear) parameter enter
%      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Swift-Hohenberg in the Hamiltonian form
lam_start=1+1*1i; % near which point do we want expansion
p=1;
A(:,:,2)=zeros(4,4);A(4,1,2)=-1; % how does the (linear) parameter enter
A(:,:,1)=[0   1  0   0;...
   -1  0  1   0;...
   0   0  0   1;...
   0   0  -1  0];...
n=4; % size of matrices
k=2; % dimension of unstable eigenspace

%  

%%%%% Cahn-Hilliard, unstable PDR at 1/4 ; spreading speed is roughly 1.62208 
%  p=1;
%  s_lin=2/3/sqrt(6)*(2+sqrt(7))*sqrt(sqrt(7)-1);
%  lam_lin=1i*(3+sqrt(7))*sqrt((2+sqrt(7))/96);
%  lam_start=1+1i; % near which point do we want expansion
%  lam_start=.5+1*1i;
%  A(:,:,2)=zeros(4,4);A(4,1,2)=-1; % how does the (linear) parameter enter
%  A(:,:,1)=[0   1  0   0;...
%     0   0  1   0;...
%     0   0  0   1;...
%     0   s_lin -1   0];...
%  n=4; % size of matrices
%  k=2; % dimension of unstable eigenspace

%%%% CPW example
%  lam_start=3;
%  A(:,:,2)=[1 0;0 -1];
%  A(:,:,1)=[0 1e2;0 0];
%  n=2;
%  k=1;


%%%% Drift-Diffusion uxx+c ux+ mu u = lam u
%  lam_start=3.1; % near which point do we want expansion
%  A(:,:,2)=zeros(2,2);A(2,1,2)=1; % how does the (linear) parameter enter
%  mu=2;c=1; % double root at mu-c^2/4
%  A(:,:,1)=[0 1; -mu -c];
%  n=2; % size of matrices
%  k=1; % dimension of unstable eigenspace


%%%% coupled Drift-Diffusion -- see mathematica notebook for values at d>2, lam=0 and nu=-1 for d<2
%% u_xx + c u_x + u +v=lam u
%% d v_xx + c v_x = lam v
%  c=2;d=3.5;
%  lam_start=4.5 + .100*1i; % near which point do we want expansion
%  A(:,:,2)=zeros(4,4);A(2,1,2)=1;A(4,3,2)=1/d; % how does the (linear) parameter enter
%  A(:,:,1)=[0 1 0 0; ...
%      -1 -c 1 0;...
%      0 0 0 1;...
%      0 0 0 -c/d];
%  n=4; % size of matrices
%  k=2; % dimension of unstable eigenspace




%%%% Laplacian in strip
%% u_xx + u_yy + c u_x + u =lam u
%  p=1;
%  lam_start=1+1i;
%  N=100;
%  c=2;
%  mu=1;
%  L=pi;
%  dx=L/N;
%  e=ones(N,1);
%  D2=spdiags([e -2*e e], -1:1,N,N);
%  %  D2(N,1)=1;D2(1,N)=1; % if periodic bc
%  D2=full(D2/dx^2);
%  A0=zeros(2*N,2*N);
%  A0(1:N,N+1:2*N)=speye(N,N);
%  A0(N+1:2*N,1:N)=-D2-mu*speye(N,N);
%  A0(N+1:2*N,N+1:2*N)=-c*speye(N,N);
%  A(:,:,1)=A0;
%  A(:,:,2)=zeros(2*N,2*N);
%  A(N+1:2*N,1:N,2)=speye(N,N);
%  n=2*N;
%  k=N;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of example definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
