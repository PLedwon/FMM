function [kzvals,xvec,kxvec,HyCoeffs] = TEsinglelayer(NumPts,period,omega,kappa,epsilonstruct,mustruct,kztolfactor)
    Delta = period/NumPts;
    xvec  = (0:(NumPts-1))*Delta;
    kxvec = 2*pi*((-NumPts/2):1:(NumPts/2-1))/period;
    

    epsXX = epsilonstruct.xx(xvec);
%     epsYY = epsilonstruct.yy(xvec);
    epsZZ = epsilonstruct.zz(xvec);
    
%     muXX = mustruct.xx(xvec);
    muYY = mustruct.yy(xvec);
%     muZZ = mustruct.zz(xvec);
    
    PDmat = 1i*diag(kxvec + kappa);
    epsXXft = fft(epsXX)/NumPts;
    muYYft = fft(muYY)/NumPts;
    invepsZZft = fft(1./epsZZ)/NumPts;
    
    epsXXmat = zeros(NumPts);
    muYYmat = zeros(NumPts);
    invepsZZmat = zeros(NumPts);
    for i=1:NumPts; 
        for j=1:NumPts;
            epsXXmat(i,j) = epsXXft(mod(i-j,NumPts)+1);
            muYYmat(i,j) = muYYft(mod(i-j,NumPts)+1);
            invepsZZmat(i,j) = invepsZZft(mod(i-j,NumPts)+1);
        end; 
    end;
    
    [Vhycoeffs,Deigvals]=eig(epsXXmat*(muYYmat*omega*omega + PDmat*invepsZZmat*PDmat));
    [kzvals,I] = sort(sqrt(diag(Deigvals)),'ascend');
    kztol = max(abs(kzvals))*kztolfactor;
    kzvals = kzfilter(kzvals,kztol);
    
    HyCoeffs = Vhycoeffs(:,I);
    
    
    Hymodestmp = zeros(size(Vhycoeffs));
    for i=1:NumPts;
        Hymodestmp(:,i) = ifft((Vhycoeffs(:,i)));
    end;
    Hymodes = Hymodestmp(:,I);
end