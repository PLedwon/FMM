function [ TransferI_Ey,TransferI_Hy ] = calculate_interface_matrix( layer_data1, layer_data2 )
%CALCULATE_INTERFACE_MATRIX Summary of this function goes here
%   Detailed explanation goes here

%Ey first

muXX1 = layer_data1.mustruct.xx(layer_data1.xvec);
muXX2 = layer_data2.mustruct.xx(layer_data2.xvec);

EyExpansionMat1 = zeros(layer_data1.NumPts);
EyExpansionMat2 = zeros(layer_data2.NumPts);

for i=1:layer_data1.NumPts;
    EyExpansionMat1(:,i) = ifft(fftshift(layer_data1.EyCoeffs(:,i)))*layer_data1.NumPts;
    EyExpansionMat2(:,i) = ifft(fftshift(layer_data2.EyCoeffs(:,i)))*layer_data2.NumPts;
end;

L1 = zeros(2*layer_data1.NumPts);
L2 = zeros(2*layer_data2.NumPts);

L1(1:layer_data1.NumPts,1:layer_data1.NumPts) = EyExpansionMat1;
L2(1:layer_data2.NumPts,1:layer_data2.NumPts) = EyExpansionMat2;

L1(1:layer_data1.NumPts,(layer_data1.NumPts+1):(2*layer_data1.NumPts)) = EyExpansionMat1;
L2(1:layer_data2.NumPts,(layer_data2.NumPts+1):(2*layer_data2.NumPts)) = EyExpansionMat2;

for i=1:layer_data1.NumPts;
    L1((layer_data1.NumPts+1):(2*layer_data1.NumPts),i) = layer_data1.kzvalsEy(i)./layer_data1.omega ./muXX1.'.*EyExpansionMat1(:,i);
    L2((layer_data2.NumPts+1):(2*layer_data2.NumPts),i) = layer_data2.kzvalsEy(i)./layer_data2.omega ./muXX2.'.*EyExpansionMat2(:,i);
    L1((layer_data1.NumPts+1):(2*layer_data1.NumPts),layer_data1.NumPts+i) = -layer_data1.kzvalsEy(i)./layer_data1.omega ./muXX1.'.*EyExpansionMat1(:,i);
    L2((layer_data2.NumPts+1):(2*layer_data2.NumPts),layer_data2.NumPts+i) = -layer_data2.kzvalsEy(i)./layer_data2.omega ./muXX2.'.*EyExpansionMat2(:,i);
end;    

clear EyExpansionMat1 EyExpansionMat2;

TransferI_Ey = L2\L1;

clear L1 L2;


%Hy second

epsXX1 = layer_data1.epsstruct.xx(layer_data1.xvec);
epsXX2 = layer_data2.epsstruct.xx(layer_data2.xvec);

HyExpansionMat1 = zeros(layer_data1.NumPts);
HyExpansionMat2 = zeros(layer_data2.NumPts);


for i=1:layer_data1.NumPts;
    HyExpansionMat1(:,i) = ifft(fftshift(layer_data1.HyCoeffs(:,i)))*layer_data1.NumPts;
    HyExpansionMat2(:,i) = ifft(fftshift(layer_data2.HyCoeffs(:,i)))*layer_data2.NumPts;
end;

L1 = zeros(2*layer_data1.NumPts);
L2 = zeros(2*layer_data2.NumPts);

L1(1:layer_data1.NumPts,1:layer_data1.NumPts) = HyExpansionMat1;
L2(1:layer_data2.NumPts,1:layer_data2.NumPts) = HyExpansionMat2;

L1(1:layer_data1.NumPts,(layer_data1.NumPts+1):(2*layer_data1.NumPts)) = HyExpansionMat1;
L2(1:layer_data2.NumPts,(layer_data2.NumPts+1):(2*layer_data2.NumPts)) = HyExpansionMat2;

for i=1:layer_data1.NumPts;
    L1((layer_data1.NumPts+1):(2*layer_data1.NumPts),i) = layer_data1.kzvalsHy(i)./layer_data1.omega ./epsXX1.'.*HyExpansionMat1(:,i);
    L2((layer_data2.NumPts+1):(2*layer_data2.NumPts),i) = layer_data2.kzvalsHy(i)./layer_data2.omega ./epsXX2.'.*HyExpansionMat2(:,i);
    L1((layer_data1.NumPts+1):(2*layer_data1.NumPts),layer_data1.NumPts+i) = -layer_data1.kzvalsHy(i)./layer_data1.omega ./epsXX1.'.*HyExpansionMat1(:,i);
    L2((layer_data2.NumPts+1):(2*layer_data2.NumPts),layer_data2.NumPts+i) = -layer_data2.kzvalsHy(i)./layer_data2.omega ./epsXX2.'.*HyExpansionMat2(:,i);
end;    

clear HyExpansionMat1 HyExpansionMat2;

TransferI_Hy = L2\L1;

clear L1 L2;

end

