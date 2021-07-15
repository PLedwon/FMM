function [ ampl_vecEy, ampl_vecHy  ] = calculate_up_input_amplitudes( layer_data, infuncEy, infuncHy, z0 )
%calculation of the input amplitudes for both given TE and TM input
%functions

Eyinc=infuncEy(layer_data.xvec);
Hyinc=infuncHy(layer_data.xvec);

EyExpansionMat = zeros(layer_data.NumPts);
HyExpansionMat = zeros(layer_data.NumPts);

for i=1:layer_data.NumPts;
        %missing kappa
        EyExpansionMat(i,:) = exp(1i*layer_data.kzvalsEy(i)*z0).*exp(1i*layer_data.kappa.*layer_data.xvec.').*ifft(fftshift(layer_data.EyCoeffs(:,i)))*layer_data.NumPts;
        HyExpansionMat(i,:) = exp(1i*layer_data.kzvalsHy(i)*z0).*exp(1i*layer_data.kappa.*layer_data.xvec.').*ifft(fftshift(layer_data.HyCoeffs(:,i)))*layer_data.NumPts;
end;

ampl_vecEy = (Eyinc/EyExpansionMat).';
ampl_vecHy = (Hyinc/HyExpansionMat).';


end

