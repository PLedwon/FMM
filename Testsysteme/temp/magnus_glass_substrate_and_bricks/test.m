%figure(1)
%imagesc(wavelength,incAngle,log10(abs(Ns8.RHy(:,:,end).'+Ns8.THy(:,:,end).'-1)))
%colorbar()

figure(2)
imagesc(wavelength,incAngle,MagnusRHy(:,:,end).')
colorbar()