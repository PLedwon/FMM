load('glass_grating_fein.mat'); close all

figure(1)
imagesc(omega,incAngle,(TEy).')
colorbar()

% figure(2)
% surf(1.5+real(ey))
% hold on
% imagesc(real(ey))
% 
