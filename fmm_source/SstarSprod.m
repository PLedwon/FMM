function [ Stot ] = SstarSprod( S1,S2 )
%SSTARSPROD Summary of this function goes here
%   Detailed explanation goes here

Stot = zeros(size(S1));

N = size(S1,1)/2;
u=1:N;
d=(N+1):(2*N);

StmpadInv = eye(N)-S1(u,d)*S2(d,u);

Stot(u,u) = S2(u,u)*(StmpadInv\S1(u,u));
Stot(u,d) = S2(u,u)*(StmpadInv\S1(u,d))*S2(d,d) + S2(u,d);
Stot(d,u) = S1(d,u) + S1(d,d)*S2(d,u)*(StmpadInv\S1(u,u));
Stot(d,d) = S1(d,d)*(S2(d,u)*(StmpadInv\S1(u,d))+eye(N))*S2(d,d);

end

