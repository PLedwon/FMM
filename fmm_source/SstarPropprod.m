function [ Stot ] = SstarPropprod( S1, kzvals,z )
%SSTARPROPPROD Summary of this function goes here
%   Detailed explanation goes here
Stot = zeros(size(S1));

N = size(S1,1)/2;

if N ~= length(kzvals)
    error('SstarPropprod input dimensions do not fit!');
end;
u=1:N;
d=(N+1):(2*N);

S2sub = diag(exp(1i*kzvals*z));

Stot(u,u) = S2sub * S1(u,u);
Stot(u,d) = S2sub*S1(u,d)*S2sub;
Stot(d,u) = S1(d,u);
Stot(d,d) = S1(d,d)*S2sub;

end

