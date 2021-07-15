function [ layer_data ] = process_single_layer( layerstack, LayerNumber )
%PROCESS_SINGLE_LAYER 

layer_data.NumPts = layerstack.NumPts;
layer_data.kappa = layerstack.kappa;
layer_data.omega = layerstack.omega;
layer_data.period = layerstack.period ;
layer_data.kztolfactor = layerstack.kztolfactor;

layer_data.thickness = layerstack.layer_thicknesses(LayerNumber);

layer_data.epsstruct = layerstack.epsstruct(LayerNumber);
layer_data.mustruct = layerstack.mustruct(LayerNumber);

%% calculate modes
[layer_data.kzvalsHy,layer_data.xvec,layer_data.kxvec,layer_data.HyCoeffs] = TEsinglelayer(layer_data.NumPts,layer_data.period,layer_data.omega,layer_data.kappa,layer_data.epsstruct,layer_data.mustruct,layer_data.kztolfactor);

[layer_data.kzvalsEy,layer_data.xvec,layer_data.kxvec,layer_data.EyCoeffs] = TMsinglelayer(layer_data.NumPts,layer_data.period,layer_data.omega,layer_data.kappa,layer_data.epsstruct,layer_data.mustruct,layer_data.kztolfactor);

end

