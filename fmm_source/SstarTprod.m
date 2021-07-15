function [ Stot ] = SstarTprod( S,T )
%SSTARTPROD Summary of this function goes here
%   Detailed explanation goes here

Stot = zeros(size(S));

N = size(S,1)/2;
u=1:N;
d=(N+1):(2*N);

Tmpadinv = T(d,u)*S(u,d)+T(d,d);
Stot(u,d) = (T(u,u)*S(u,d) + T(u,d))/Tmpadinv;
Stot(u,u) = T(u,u)*S(u,u) - Stot(u,d)*T(d,u)*S(u,u);
Stot(d,d) = S(d,d)/Tmpadinv;
Stot(d,u) = S(d,u) -Stot(d,d)*T(d,u)*S(u,u);

end

