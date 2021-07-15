function [kzvals,xvec,kxvec,EyCoeffs] = TMsinglelayer(NumPts,period,omega,kappa,epsilonstruct,mustruct,kztolfactor)

    Delta = period/NumPts;
    xvec  = (0:(NumPts-1))*Delta;
    kxvec = 2*pi*((-NumPts/2):1:(NumPts/2-1))/period;
    

%     epsXX = epsilonstruct.xx(xvec);
    epsYY = epsilonstruct.yy(xvec);
%     epsZZ = epsilonstruct.zz(xvec);
    
    muXX = mustruct.xx(xvec);
%     muYY = mustruct.yy(xvec);
    muZZ = mustruct.zz(xvec);

    PDmat = 1i*diag(kxvec + kappa);
    epsYYft = fft(epsYY)/NumPts;
    muXXft = fft(muXX)/NumPts;
    invmuZZft = fft(1./muZZ)/NumPts;
    
    epsYYmat = zeros(NumPts);
    muXXmat = zeros(NumPts);
    invmuZZmat = zeros(NumPts);
    for i=1:NumPts; 
        for j=1:NumPts;
            epsYYmat(i,j) = epsYYft(mod(i-j,NumPts)+1);
            muXXmat(i,j) = muXXft(mod(i-j,NumPts)+1);
            invmuZZmat(i,j) = invmuZZft(mod(i-j,NumPts)+1);
        end; 
    end;
    
    [Veycoeffs,Deigvals]=eig(muXXmat*(epsYYmat*omega*omega + PDmat*invmuZZmat*PDmat));
    [kzvals,I] = sort(sqrt(diag(Deigvals)));
    kztol = max(abs(kzvals))*kztolfactor;
    kzvals = kzfilter(kzvals,kztol);
    
    EyCoeffs = Veycoeffs(:,I);
end
