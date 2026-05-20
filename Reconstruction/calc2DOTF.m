function T2D = calc2DOTF(lambda,NA,FOVx,pxN,P,S)
dx = FOVx/pxN; % Image pixel size in µm
kx_max = 0.5/dx; % Max spatial frequency in lines/µm after fft of image
dkx = 2*kx_max/pxN; % Frequency step per pixel in lines/µm of image fft
fc_NA = 2*NA/lambda; % Cut-off frequency from MO in lines/µm
mask_fc = fc_NA/dkx;
kxBins = linspace(-kx_max,kx_max-dkx,pxN); % image frequency space k axis

T2D = zeros(pxN,pxN);

for i = 1:pxN/2
    tic
    for j = 1:pxN/2
        qx = kxBins(i); % Select a qx
        qy = kxBins(j); % Select a qy
        % Only continue if shifted pupils would overlap as the contribution
        % is zero otherwise
        if(sqrt((qx/2)^2 + (qy/2)^2) <= (ceil(mask_fc))*dkx)
            P_plus  = imtranslate(P,  [qx/2/dkx, qy/2/dkx]); % Shift pupil
            P_minus = rot90(P_plus,2); % rotate 180° due to symmetry
            S_plus  = imtranslate(S,  [qx/2/dkx, qy/2/dkx]); % Shift illumination distribution
            S_minus = -rot90(S_plus,2); % rotate 180° and negate due to symmetry
            SP = (S_plus-S_minus).*P_plus.*P_minus; % Integrand without delta function (T2D)
            T2D(i,j) = sum(SP,'all'); % Integrate SP to calculate T2D at (qx,qz)            
        end
    end
    toc
    fprintf("Done with sub-integration %d of %d!\n",i,pxN/2);
end
% Calculate the rest of T2D via mirroring according to its symmetry
T2D = T2D + fliplr(T2D);
T2D = T2D - flipud(T2D);
end