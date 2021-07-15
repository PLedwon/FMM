function [ P0incEy,P0scatEy,PendincEy,PendscatEy,P0incHy,P0scatHy,PendincHy,PendscatHy,Ex,Ey,Ez,Hx,Hy,Hz,diffraction_order_data] = smat_with_modes( layerstack,xvec,zvec )
%SMAT_WITH_MODES Summary of this function goes here
%   Detailed explanation goes here

%% process all layers
if length(layerstack.layer_thicknesses) ~= length(layerstack.epsstruct) || length(layerstack.layer_thicknesses) ~= length(layerstack.mustruct);
    error('layer_thicknessvector and eps or mu structure are not the same length!');
end;

num_of_layers = length(layerstack.layer_thicknesses);

fprintf('prepare layers: \n');
tic;
%layerstack.kappa
for i=1:num_of_layers;
    layer_data(i) = process_single_layer( layerstack, i );
end;
toc;

%% determine initial amplitudes
fprintf('get initial amplitudes: \n');
tic;
[ u1Ey, u1Hy  ] = calculate_up_input_amplitudes( layer_data(1), layerstack.UinfuncEy, layerstack.UinfuncHy, 0 );
[ dNEy, dNHy  ] = calculate_down_input_amplitudes( layer_data(num_of_layers), layerstack.DinfuncEy, layerstack.DinfuncHy, 0 );
toc;
%% initialize amplitude structs
uampls(1).Ey = u1Ey;
uampls(1).Hy = u1Hy;

dampls(num_of_layers).Ey = dNEy;
dampls(num_of_layers).Hy = dNHy;

%% quasi-recursive algorithmn over all layers (not to save all sub-S-Matrix results)

upper_uampls1_Ey = exp(1i*layer_data(1).kzvalsEy*layer_data(1).thickness).*uampls(1).Ey;
upper_uampls1_Hy = exp(1i*layer_data(1).kzvalsHy*layer_data(1).thickness).*uampls(1).Hy;

N = layerstack.NumPts;

fprintf('determine modes:\n');
tic;
for j = num_of_layers:-1:2;
    %% propagate to down ampls to 'lower' layer bdry or determine them from previous loop
    if j == num_of_layers;
%         lower_damplsj_Ey = exp(1i*layer_data(j).kzvalsEy*layer_data(j).thickness).*dNEy;
% exp(1i*layer_data(j).kzvalsHy*layer_data(j).thickness).*
        dampls(j).Ey = dNEy;
        dampls(j).Hy = dNHy;
    elseif j < num_of_layers;
        [ invTransferI_Ey,invTransferI_Hy ] = calculate_interface_matrix( layer_data(j+1), layer_data(j) );
        dampls(j).Ey = invTransferI_Ey(1:N,(N+1):(2*N))*uampls(j+1).Ey  + invTransferI_Ey((N+1):(2*N),(N+1):(2*N))*(exp(1i*layer_data(j+1).kzvalsEy*layer_data(j+1).thickness).*dampls(j+1).Ey);
        dampls(j).Hy = invTransferI_Hy(1:N,(N+1):(2*N))*uampls(j+1).Hy  + invTransferI_Hy((N+1):(2*N),(N+1):(2*N))*(exp(1i*layer_data(j+1).kzvalsHy*layer_data(j+1).thickness).*dampls(j+1).Hy);
    else
        error('index exceeds number of layers!');
    end;
    
    %% calculate S-Matrix for specific setup
    StotEy = eye(2*layerstack.NumPts);
    StotHy = eye(2*layerstack.NumPts);
    for i=1:j-2;
        [ TransferI_Ey,TransferI_Hy ] = calculate_interface_matrix( layer_data(i), layer_data(i+1) );
        StotEy = SstarTprod(StotEy,TransferI_Ey);
        StotHy = SstarTprod(StotHy,TransferI_Hy);
        StotEy = SstarPropprod(StotEy, layer_data(i+1).kzvalsEy , layer_data(i+1).thickness); 
        StotHy = SstarPropprod(StotHy, layer_data(i+1).kzvalsHy , layer_data(i+1).thickness);
    end;
    [ TransferI_Ey,TransferI_Hy ] = calculate_interface_matrix( layer_data(j-1), layer_data(j) );
    StotEy = SstarTprod(StotEy,TransferI_Ey);
    StotHy = SstarTprod(StotHy,TransferI_Hy);
    
    %% calculate missing up and down ampls for the missing layer
    
    uampls(j).Ey = StotEy(1:N,1:N)*upper_uampls1_Ey + StotEy(1:N,(N+1):(2*N))*(exp(1i*layer_data(j).kzvalsEy*layer_data(j).thickness).*dampls(j).Ey);
    uampls(j).Hy = StotHy(1:N,1:N)*upper_uampls1_Hy+ StotHy(1:N,(N+1):(2*N))*(exp(1i*layer_data(j).kzvalsHy*layer_data(j).thickness).*dampls(j).Hy);
    
    dampls(1).Ey = StotEy((N+1):(2*N),1:N)*upper_uampls1_Ey + StotEy((N+1):(2*N),(N+1):(2*N))*(exp(1i*layer_data(j).kzvalsEy*layer_data(j).thickness).*dampls(j).Ey);
    dampls(1).Hy = StotHy((N+1):(2*N),1:N)*upper_uampls1_Hy + StotHy((N+1):(2*N),(N+1):(2*N))*(exp(1i*layer_data(j).kzvalsHy*layer_data(j).thickness).*dampls(j).Hy);
        
end;
clear StotEy StotHy;

%% determine fields for each layer
Ex = zeros(length(xvec),length(zvec));
Ey = zeros(length(xvec),length(zvec));
Ez = zeros(length(xvec),length(zvec));
Hx = zeros(length(xvec),length(zvec));
Hy = zeros(length(xvec),length(zvec));
Hz = zeros(length(xvec),length(zvec));
% lower half plane

[ Extmp,Eytmp,Eztmp,Hxtmp,Hytmp,Hztmp ] = get_layer_modes( layer_data(1) , xvec );

lowerzbdry = 0;
upperzbdry = layerstack.layer_thicknesses(1);

for k = 1:length(zvec);
    if zvec(k) <= upperzbdry
        % Ey related fields
        current_Ey_uampl = exp(1i*layer_data(1).kzvalsEy*(zvec(k)-lowerzbdry)).*uampls(1).Ey;
        current_Ey_dampl = exp(1i*layer_data(1).kzvalsEy*(upperzbdry-zvec(k))).*dampls(1).Ey;
        
        Hx(:,k) = Hxtmp*current_Ey_uampl - Hxtmp*current_Ey_dampl;
        Ey(:,k) = Eytmp*current_Ey_uampl + Eytmp*current_Ey_dampl;
        Hz(:,k) = Hztmp*current_Ey_uampl + Hztmp*current_Ey_dampl;
        
        % Hy related fields
        current_Hy_uampl = exp(1i*layer_data(1).kzvalsHy*(zvec(k)-lowerzbdry)).*uampls(1).Hy;
        current_Hy_dampl = exp(1i*layer_data(1).kzvalsHy*(upperzbdry-zvec(k))).*dampls(1).Hy;
        
        Ex(:,k) = Extmp*current_Hy_uampl - Extmp*current_Hy_dampl;
        Hy(:,k) = Hytmp*current_Hy_uampl + Hytmp*current_Hy_dampl;
        Ez(:,k) = Eztmp*current_Hy_uampl + Eztmp*current_Hy_dampl;
    end;
end;

% mid planes

for j = 2:num_of_layers-1;
    [ Extmp,Eytmp,Eztmp,Hxtmp,Hytmp,Hztmp ] = get_layer_modes( layer_data(j) , xvec );
    
    lowerzbdry = sum(layerstack.layer_thicknesses(1:(j-1)));
    upperzbdry = sum(layerstack.layer_thicknesses(1:j));
    
    for k = 1:length(zvec);
        if zvec(k) > lowerzbdry && zvec(k) <= upperzbdry;
            % Ey related fields
            current_Ey_uampl = exp(1i*layer_data(j).kzvalsEy*(zvec(k)-lowerzbdry)).*uampls(j).Ey;
            current_Ey_dampl = exp(1i*layer_data(j).kzvalsEy*(upperzbdry - zvec(k))).*dampls(j).Ey;

            Hx(:,k) = Hxtmp*current_Ey_uampl - Hxtmp*current_Ey_dampl;
            Ey(:,k) = Eytmp*current_Ey_uampl + Eytmp*current_Ey_dampl;
            Hz(:,k) = Hztmp*current_Ey_uampl + Hztmp*current_Ey_dampl;

            % Hy related fields
            current_Hy_uampl = exp(1i*layer_data(j).kzvalsHy*(zvec(k)-lowerzbdry)).*uampls(j).Hy;
            current_Hy_dampl = exp(1i*layer_data(j).kzvalsHy*(upperzbdry - zvec(k))).*dampls(j).Hy;

            Ex(:,k) = Extmp*current_Hy_uampl - Extmp*current_Hy_dampl;
            Hy(:,k) = Hytmp*current_Hy_uampl + Hytmp*current_Hy_dampl;
            Ez(:,k) = Eztmp*current_Hy_uampl + Eztmp*current_Hy_dampl;
        end;
    end;
end;
% upper half plane
[ Extmp,Eytmp,Eztmp,Hxtmp,Hytmp,Hztmp ] = get_layer_modes( layer_data(num_of_layers) , xvec );
lowerzbdry = sum(layerstack.layer_thicknesses(1:(num_of_layers-1)));
upperzbdry = sum(layerstack.layer_thicknesses(1:num_of_layers));
for k = 1:length(zvec);
    if zvec(k) > lowerzbdry
        % Ey related fields
        current_Ey_uampl = exp(1i*layer_data(num_of_layers).kzvalsEy*(zvec(k)-lowerzbdry)).*uampls(num_of_layers).Ey;
        current_Ey_dampl = exp(1i*layer_data(num_of_layers).kzvalsEy*(upperzbdry - zvec(k))).*dampls(num_of_layers).Ey;

        Hx(:,k) = Hxtmp*current_Ey_uampl - Hxtmp*current_Ey_dampl;
        Ey(:,k) = Eytmp*current_Ey_uampl + Eytmp*current_Ey_dampl;
        Hz(:,k) = Hztmp*current_Ey_uampl + Hztmp*current_Ey_dampl;

        % Hy related fields
        current_Hy_uampl = exp(1i*layer_data(num_of_layers).kzvalsHy*(zvec(k)-lowerzbdry)).*uampls(num_of_layers).Hy;
        current_Hy_dampl = exp(1i*layer_data(num_of_layers).kzvalsHy*(upperzbdry - zvec(k))).*dampls(num_of_layers).Hy;

        Ex(:,k) = Extmp*current_Hy_uampl - Extmp*current_Hy_dampl;
        Hy(:,k) = Hytmp*current_Hy_uampl + Hytmp*current_Hy_dampl;
        Ez(:,k) = Eztmp*current_Hy_uampl + Eztmp*current_Hy_dampl;
    end;
end;
toc;
%% for power transport calculation
%% calculate S-Matrix again (up to desired thickness in upper halfspace, i.e. avoid last T-Matrix multiplication)
fprintf('do final S-Matrix:\n');
tic;
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
toc;
fprintf('\n');
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
    diffraction_order_data.outAngle0(i) = real(180./pi*asin((layer_data(1).kxvec(i)+layer_data(1).kappa)./layer_data(1).omega./sqrt(layer_data(1).epsstruct.xx(0))));% ...
        %- imag(180./pi*asin((layer_data(1).kxvec(i)+layer_data(1).kappa)./layer_data(1).omega./sqrt(layer_data(1).epsstruct.xx(0))));
    diffraction_order_data.outAngleend(i) = real(180./pi*asin((layer_data(num_of_layers).kxvec(i)+layer_data(num_of_layers).kappa)./layer_data(num_of_layers).omega./sqrt(layer_data(num_of_layers).epsstruct.xx(0))));% ...
        %- imag(180./pi*asin((layer_data(num_of_layers).kxvec(i)+layer_data(num_of_layers).kappa)./layer_data(num_of_layers).omega./sqrt(layer_data(num_of_layers).epsstruct.xx(0))));
end;


end


