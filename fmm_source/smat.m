function [ P0incEy,P0scatEy,PendincEy,PendscatEy,P0incHy,P0scatHy,PendincHy,PendscatHy,diffraction_order_data,StotEy,StotHy] = smat( layerstack )

%% process all layers
if length(layerstack.layer_thicknesses) ~= length(layerstack.epsstruct) || length(layerstack.layer_thicknesses) ~= length(layerstack.mustruct);
    error('layer_thicknessvector and eps or mu structure are not the same length!');
end;

num_of_layers = length(layerstack.layer_thicknesses);

% fprintf('prepare layers: \n');
% tic;
for i=1:num_of_layers;
    layer_data(i) = process_single_layer( layerstack, i );
end;
% toc;

%% determine initial amplitudes
% fprintf('get initial amplitudes: \n');
% tic;
[ u1Ey, u1Hy  ] = calculate_up_input_amplitudes( layer_data(1), layerstack.UinfuncEy, layerstack.UinfuncHy, 0 );
[ dNEy, dNHy  ] = calculate_down_input_amplitudes( layer_data(num_of_layers), layerstack.DinfuncEy, layerstack.DinfuncHy, 0 );
% toc;

%% calculate S-Matricies
% fprintf('do S-Matrix:\n');
% tic;
StotEy = eye(2*layerstack.NumPts);
StotHy = eye(2*layerstack.NumPts);
StotEy = SstarPropprod( StotEy , layer_data(1).kzvalsEy , layer_data(1).thickness);
StotHy = SstarPropprod( StotHy , layer_data(1).kzvalsHy , layer_data(1).thickness);
for i=2:num_of_layers;
    [ TransferI_Ey,TransferI_Hy ] = calculate_interface_matrix( layer_data(i-1), layer_data(i) );
    StotEy = SstarTprod(StotEy,TransferI_Ey);
    StotHy = SstarTprod(StotHy,TransferI_Hy);
    StotEy = SstarPropprod(StotEy, layer_data(i).kzvalsEy , layer_data(i).thickness); 
    StotHy = SstarPropprod(StotHy, layer_data(i).kzvalsHy , layer_data(i).thickness);
end;
% toc;
% fprintf('\n');
%% calculation of the amplitude vectors
N = layerstack.NumPts;

uNEy = StotEy(1:N,1:N)*u1Ey + StotEy(1:N,(N+1):(2*N))*dNEy;
d1Ey = StotEy((N+1):(2*N),1:N)*u1Ey + StotEy((N+1):(2*N),(N+1):(2*N))*dNEy;

uNHy = StotHy(1:N,1:N)*u1Hy + StotHy(1:N,(N+1):(2*N))*dNHy;
d1Hy = StotHy((N+1):(2*N),1:N)*u1Hy + StotHy((N+1):(2*N),(N+1):(2*N))*dNHy;

%% calculation of the corresponding power transported
% assumption of of non structered materials in 1st and last layer, i.e.
% const in x-direction.

%field expansion coefficients of the regions the poyntingvector shall
%beevaluated
EyCoeffs1 = layer_data(1).EyCoeffs;
HxCoeffs1 = zeros(size(EyCoeffs1));
for m = 1:size(EyCoeffs1,2);
    HxCoeffs1(:,m) = - layer_data(1).kzvalsEy(m)./layer_data(1).omega./layer_data(1).mustruct.xx(0) * EyCoeffs1(:,m);
end;

EyCoeffsN = layer_data(num_of_layers).EyCoeffs;
HxCoeffsN = zeros(size(EyCoeffsN));
for m = 1:size(EyCoeffsN,2);
    HxCoeffsN(:,m) = - layer_data(num_of_layers).kzvalsEy(m)./layer_data(num_of_layers).omega./layer_data(num_of_layers).mustruct.xx(0) * EyCoeffsN(:,m);
end;


HyCoeffs1 = layer_data(1).HyCoeffs;
ExCoeffs1 = zeros(size(EyCoeffs1));
for m = 1:size(HyCoeffs1,2);
    ExCoeffs1(:,m) = layer_data(1).kzvalsHy(m)./layer_data(1).omega./layer_data(1).epsstruct.xx(0) * HyCoeffs1(:,m);
end;


HyCoeffsN = layer_data(num_of_layers).HyCoeffs;
ExCoeffsN = zeros(size(EyCoeffsN));
for m = 1:size(HyCoeffsN,2);
    ExCoeffsN(:,m) = layer_data(num_of_layers).kzvalsHy(m)./layer_data(num_of_layers).omega./layer_data(num_of_layers).epsstruct.xx(0) * HyCoeffsN(:,m);
end;

% from incedent fields

P0incEy = - 0.5*real(u1Ey'*HxCoeffs1'*EyCoeffs1*u1Ey); 
P0incHy = 0.5*real(u1Hy'*ExCoeffs1'*HyCoeffs1*u1Hy);
PendincEy = 0.5*real(dNEy'*HxCoeffsN'*EyCoeffsN*dNEy);
PendincHy = - 0.5*real(dNHy'*ExCoeffsN'*HyCoeffsN*dNHy);

% from scattered fields

P0scatEy = 0.5*real(d1Ey'*HxCoeffs1'*EyCoeffs1*d1Ey);
P0scatHy = - 0.5*real(d1Hy'*ExCoeffs1'*HyCoeffs1*d1Hy);
PendscatEy = - 0.5*real(uNEy'*HxCoeffsN'*EyCoeffsN*uNEy);
PendscatHy = + 0.5*real(uNHy'*ExCoeffsN'*HyCoeffsN*uNHy);

% get diffraction order data:
diffraction_order_data.order_number = (-layerstack.NumPts/2):1:(layerstack.NumPts/2-1);
for i=1:layerstack.NumPts;
    diffraction_order_data.P0scatEy(i) = 0.5*real(d1Ey'*HxCoeffs1(i,:)'*EyCoeffs1(i,:)*d1Ey);
    diffraction_order_data.P0scatHy(i) = - 0.5*real(d1Hy'*ExCoeffs1(i,:)'*HyCoeffs1(i,:)*d1Hy);
    diffraction_order_data.PendscatEy(i) = - 0.5*real(uNEy'*HxCoeffsN(i,:)'*EyCoeffsN(i,:)*uNEy);
    diffraction_order_data.PendscatHy(i) = + 0.5*real(uNHy'*ExCoeffsN(i,:)'*HyCoeffsN(i,:)*uNHy);
    if abs((layer_data(1).kxvec(i)+layer_data(1).kappa)./layer_data(1).omega./sqrt(layer_data(1).epsstruct.xx(0))) > 1
        diffraction_order_data.is_evanescent0(i) = 1;
    else
        diffraction_order_data.is_evanescent0(i) = 0;
    end;
    if abs((layer_data(num_of_layers).kxvec(i)+layer_data(num_of_layers).kappa)./layer_data(num_of_layers).omega./sqrt(layer_data(num_of_layers).epsstruct.xx(0))) > 1
        diffraction_order_data.is_evanescentend(i) = 1;
    else
        diffraction_order_data.is_evanescentend(i) = 0;
    end;
    diffraction_order_data.outAngle0(i) = real(180./pi*asin((layer_data(1).kxvec(i)+layer_data(1).kappa)./layer_data(1).omega./sqrt(layer_data(1).epsstruct.xx(0))));
    diffraction_order_data.outAngleend(i) = real(180./pi*asin((layer_data(num_of_layers).kxvec(i)+layer_data(num_of_layers).kappa)./layer_data(num_of_layers).omega./sqrt(layer_data(num_of_layers).epsstruct.xx(0))));
end;


%%
% P0scat = d1Ey;
% Pendscat = uNEy;
% 
% 
% 
% P0inc = u1Ey;
% Pendinc = dNEy;
%% post processing
% x = -500:1:1500;
% x = layer_data(1).xvec;
% [ Ex,Ey,Ez,Hx,Hy,Hz ] = get_layer_modes( layer_data(1) , x);
% %%
% figure(1);
% plot(x,real(Ey*u1Ey),'.-')


end

