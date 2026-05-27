%% Load the image z-stack slice by slice and preprocess
% Modify according to your file and dataformat. 
% This script also includes an extra SHG channel. Delete if not needed. 

clear('PGX_norm','PGY_norm','PGX_data','PGY_data','SUM_data' ... 
    ,'im1save','im2save','im3save','im4save','i','samplerate' ...
    ,'PGX','PGY');
load(sprintf('parameters.mat'));
load(sprintf('Zstack_IMG%d.mat',1));

% Used to cut out the image edges in order to remove scanning artifacts
% that could disturb the reconstruction.
ROI_cutout = 450;
ogImsize = length(im1save);
if ROI_cutout < ogImsize
    imSize = ROI_cutout;
else
    imSize = ogImsize;
end
Zstack_nr = length(Zpositions);

SUM = zeros(imSize, imSize, Zstack_nr);
PGX_raw = zeros(imSize, imSize, Zstack_nr);
PGY_raw = zeros(imSize, imSize, Zstack_nr);
SHG = zeros(imSize, imSize, Zstack_nr);
for i = 1:Zstack_nr
    fprintf("Loading slice %d of %d!\n",i,Zstack_nr);
    load(sprintf('Zstack_IMG%d.mat',i));
    if ROI_cutout < ogImsize        
        PGX_raw(:,:,i) = imcrop(-im1save,[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
        PGY_raw(:,:,i) = imcrop(-im2save,[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
        SUM(:,:,i) = imcrop(abs(im3save),[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
        SHG(:,:,i) = imcrop(abs(im4save),[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
    else       
        PGX_raw(:,:,i) = -im1save;
        PGY_raw(:,:,i) = -im2save;
        SUM(:,:,i) = abs(im3save);
        SUM(:,:,i) = abs(im4save);
    end
end
    
PGX(:,:,:) = PGX_raw(:,:,:) ./ SUM(:,:,:);
PGY(:,:,:) = PGY_raw(:,:,:) ./ SUM(:,:,:);

% Optional. Subtract a heavily smoothed copy to remove possible background gradients
for i = 1:Zstack_nr   
    PGX(:,:,i) = PGX(:,:,i) - imgaussfilt(PGX(:,:,i),imSize/5);
    PGY(:,:,i) = PGY(:,:,i) - imgaussfilt(PGY(:,:,i),imSize/5);
end

fprintf('Z-stack loaded!\n');
clear('PGX_raw','PGY_raw','im1save','im2save','im3save','im4save','i','samplerate','Z0');

figure;
orthosliceViewer(SUM);
colormap gray
figure;
orthosliceViewer(PGX);
colormap gray
figure;
orthosliceViewer(PGY);
colormap gray
figure;
orthosliceViewer(SHG);
colormap gray
%% Optional rescaling for up/downsampling 
% Used to resize the lateral frequency space. This is needed if you subsample 
% your image, i.e, your resolution is below the resolution of your MOs NA.
% If this happens, the image frequency space would be smaller than the
% pupil function of your MO and will lead to errors in the following
% reconstruction. 
scaleFactor = 2;
PGX_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
PGY_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
SUM_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
SHG_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
for i = 1:Zstack_nr       
    PGX_new(:,:,i) = imresize(PGX(:,:,i),scaleFactor);
    PGY_new(:,:,i) = imresize(PGY(:,:,i),scaleFactor);
    SUM_new(:,:,i) = imresize(SUM(:,:,i),scaleFactor);
    SHG_new(:,:,i) = imresize(SHG(:,:,i),scaleFactor);
end
PGX = PGX_new;
PGY = PGY_new;
SUM = SUM_new;
SHG = SHG_new;
imSize = imSize*scaleFactor;
clear('PGX_new','PGY_new','SUM_new','SHG_new');
figure;
orthosliceViewer(PGX);
colormap gray
%% 
% Microscope and image parameters:
lambda = 1.050; % Wavelength in µm
NA = 0.8; % NA of microscope objective
FOVx = 276.93; % Field of view in the aquired images in µm. Square images are assumed!
FOVx = FOVx * ROI_cutout/ogImsize; % Correct for the ROI cutout, if used
FOVz = Zstack_size; % Axial field of view in µm
pxN = imSize; % Lateral pixel number
pzN = Zstack_nr; % Axial pixel number

% Calculate lateral frequency axis
dx = FOVx/pxN; % Lateral pixel size in µm
kx_max = 0.5/dx; % Max lateral spatial frequency in lines/µm
dkx = 2*kx_max/pxN; % Lateral frequency step per pixel in lines/µm

fc_NAp = 2*NA/lambda; % Cut-off frequency from illumination MO in lines/µm
mask_fc = fc_NAp/dkx; % Pupil mask size

% Calculate axial frequency axis
dz = FOVz/pzN; % Axial image pixel size in µm
kz_max = 0.5/dz; % Max axial spatial frequency in lines/µm
dkz = 2*kz_max/pzN; % Axial frequency step per pixel in lines/µm

kxBins = linspace(-kx_max,kx_max-dkx,pxN); % Lateral image frequency space k axis
kzBins = linspace(-kz_max,kz_max-dkz,pzN); % Axial image frequency space k axis

% Calculate pupil function
% alpha =  BFA / w, ratio between the 1/e^2 beam size w and the 
% back focal aparture diameter (BFA)
% alpha = 0 defines a gaussian light source with infinite size,
% essentially plane waves at the back focal aparture.
alpha = 0.75; 
P = calcPupil(lambda,NA,FOVx,pxN,alpha);
%P_NA1 = calcPupil(lambda,1,FOVx,pxN,0);

% Display the pupil function to check size
img = P;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap gray
colorbar

%figure;
%plot(P(:,451));
%figform;

% Calculate DPLD (analog to source function S) from the simulation file
sim_data_path = 'G:\User\Nick Lemberger\Promotion\Paper\2025 --- 4Q-Epi Detector\GitHub Upload\DPLD_ChickenMuscle.jdat';
DPLD = calcDPLD(lambda,FOVx,pxN,sim_data_path);
DPLD = imgaussfilt(DPLD, pxN/50); % Smooth to reduce MC simulation noise


img = DPLD;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap jet
colorbar
%% Calculate T2D and T3D
% Calculating the optical transferfunction. This is the most time consuming step. 
% gridRes is the grid resolution at which the delta function in the T3D
% integral is evaluated in sperical coordinates. 512 is sufficiently precise
% and fast. 
gridRes = 512;
OTF = calc2D3DOTF(lambda,NA,FOVx,FOVz,pxN,pzN,gridRes,P,DPLD);

T2D = OTF.T2D;
T3D = OTF.T3D;

% Save the calculated optical transferfunctions as they only need to be
% calculated once for a given field of view, resolution, NA, and DPLD.
save('OTF_justCalced.mat', 'T2D', 'T3D')

img = T2D;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap jet
colorbar

img = sum(T3D,3);
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap jet
colorbar

figure;
orthosliceViewer(permute(-T3D,[2,1,3]));
colormap jet
colorbar

%% Optional prefiltering. Remove stuff outside of NA
for slice = 1:Zstack_nr
   PGX_ = PGX(:,:,slice);
   PGY_ = PGY(:,:,slice);
   %SHG_ = SHG(:,:,slice);
   % Apply NA of illumination MO as mask to remove out of band noise
   F_PGX_ = fftshift(fft2(PGX_)).* (P>0);
   F_PGY_ = fftshift(fft2(PGY_)).* (P>0);
   %F_SHG_ = fftshift(fft2(SHG_)).* (P>0);
   
   % Remove some specific noise bands
   % Carefully inspect the image fft.
   % These following frequency components (for 900x900px) should be 0 
   % anyway as phase gradient images have no DC component. 
   F_PGX_(:,451) = 0; 
   F_PGY_(:,451) = 0;
   F_PGX_(451,:) = 0;
   F_PGY_(451,:) = 0;
   %F_SHG_(451,:) = 0;
   %F_SHG_(451,:) = 0;
   % Specific noise bands that were notable in our setup
   %F_PGX_(:,657) = 0;
   %F_PGY_(:,657) = 0;
   %F_PGX_(:,245) = 0;
   %F_PGY_(:,245) = 0;
   
   PGX(:,:,slice) = real(ifft2(ifftshift(F_PGX_)));
   PGY(:,:,slice) = real(ifft2(ifftshift(F_PGY_)));
   %SHG(:,:,slice) = real(ifft2(ifftshift(F_SHG_)));
end
clear('PGX_','PGY_','SHG_','F_PGX_','F_PGY_','F_SHG_');

%% Reconstruct Z-Stack with T3D
alpha =5.15e-6; % Regularization parameter for direct deconvolution
% Either tune alpha with a known sample to match to the expected 
% index values or use an automatic optimization algorithm

% n0 is the background or mean refractive index across your focal plane,
% not necessarily your scattering medium!
%n0 = 1.5037; % Background index 518F immersion oil
n0 = 1.46; % Adipose chicken tissue
%n0 = 1.37; % Chicken muscle tissue
%n0 = 1.36; % Brain

% Calculate the scattering ponential V
% The sign of the PGX, PGY and T3D depends on the polarity of your phase
% gradient signals. The order of PGX and PGY may be switched depending on
% your exact setup. 
V_T3D = applyT3D(PGX,PGY, lambda/(4*pi) * T3D, alpha);
k = 2*pi*n0 / lambda;
%n = sqrt(1-V_T3D/k^2)*n0; % Full result
n = n0-n0*V_T3D/k^2; % Linear approximation


figure;
orthosliceViewer(n);
%caxis([1.27, 1.45])
colormap gray
colorbar


%% Reconstruct a single 2D image
% Choose alpha for best image quality. Chosing a large alpha will result
% in more high frequency content and sharper images, but low frequencies
% will be damped and a noticable drop of intensity towards the middle of
% objects becomes apparent. This is a common problem when using only a 2D
% reconstruction. 
alpha =5e-3; % 
slice = 50;
phiT2D = applyT2D(-PGX(:,:,slice),-PGY(:,:,slice),-T2D,alpha);

img = -phiT2D;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap gray
colorbar
