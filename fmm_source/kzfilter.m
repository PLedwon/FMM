function [ kz_new ] = kzfilter( kz_old,kztol )
%KZFILTER to determine up wards damped waves
kz_new = zeros(size(kz_old));

for i=1:length(kz_old);
    if real(kz_old(i)) >= kztol;
       if imag(kz_old(i)) < -kztol;
           kz_new(i) = -kz_old(i);
       else
           kz_new(i) = kz_old(i);
       end;
    elseif real(kz_old(i)) <= -kztol;
       if imag(kz_old(i)) <= kztol;
           kz_new(i) = -kz_old(i);
       else
           kz_new(i) = kz_old(i);
       end;
    else
       if imag(kz_old(i)) > kztol; 
           kz_new(i) = kz_old(i);
       elseif imag(kz_old(i)) < - kztol;
           kz_new(i) = -kz_old(i);
       else
           kz_new(i) = real(kz_old(i)) + 1i*(1+eps(1))*kztol; 
       end;
    end;     
end;
return;

end

