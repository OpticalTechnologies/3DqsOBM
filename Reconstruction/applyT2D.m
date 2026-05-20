function phi = applyT2D(PGX,PGY,T2D,alpha)
FPGX = fftshift(fft2(PGX));
FPGY = fftshift(fft2(PGY));
% Wiener filtering for deconvolution with T2Dx and T2Dy
denom = (abs(T2D').^2 + abs(T2D).^2 + alpha);
Fphi = (conj(1i*T2D').*FPGX + conj(1i*T2D).*FPGY) ./ denom;
phi = real(ifft2(ifftshift(Fphi))); 
end