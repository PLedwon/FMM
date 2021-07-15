function [ Ex,Ey,Ez,Hx,Hy,Hz ] = get_layer_modes( layer_data , x )
%GET_LAYER_MODES Summary of this function goes here
%   Detailed explanation goes here

Ex = zeros(length(x),size(layer_data.kxvec,1));

Ez = zeros(length(x),size(layer_data.kxvec,1));
Hx = zeros(length(x),size(layer_data.kxvec,1));

Hz = zeros(length(x),size(layer_data.kxvec,1));

%matrix for discrete FT
DFTmatW = zeros(length(x),size(layer_data.kxvec,1));
%matrix for first derivative discrete FT
PDxDFTmatW = zeros(length(x),size(layer_data.kxvec,1));
for i=1:length(layer_data.kxvec);
    for j=1:length(x);
        DFTmatW(j,i) = exp(1i*(layer_data.kxvec(i) + layer_data.kappa)*x(j));
        PDxDFTmatW(j,i) = 1i*(layer_data.kxvec(i) + layer_data.kappa)*DFTmatW(j,i);
    end;
end;

Ey = DFTmatW*layer_data.EyCoeffs;
Hy = DFTmatW*layer_data.HyCoeffs;

for i=1:size(layer_data.EyCoeffs,2);
    Ex(:,i) = layer_data.kzvalsHy(i)./layer_data.omega./layer_data.epsstruct.xx(x).'.*(DFTmatW*layer_data.HyCoeffs(:,i));
    Hx(:,i) = -layer_data.kzvalsEy(i)./layer_data.omega./layer_data.mustruct.xx(x).'.*(DFTmatW*layer_data.EyCoeffs(:,i));
    
    Ez(:,i) = 1i./layer_data.omega./layer_data.epsstruct.zz(x).'.*(PDxDFTmatW*layer_data.HyCoeffs(:,i));
    Hz(:,i) = -1i./layer_data.omega./layer_data.mustruct.zz(x).'.*(PDxDFTmatW*layer_data.EyCoeffs(:,i));
end;





end

