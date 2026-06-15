% This is a simple Matlab script to simulate the detected-photon launch
% distribution (DPLD) required for phase or index retireval in qsOBM with
% the S11142-10 photodetector.
% Author: Nick Sidney Lemberger, 19.05.2026
%% Add the path to your MXCLab installation
addpath(genpath('G:\User\Nick Lemberger\MCXStudio'));
%savepath
%%
mcxlab('gpuinfo')

%% General description of the simulation. Scaling: voxel unit size = 0.04mm
% MCX Simulation for the S11142-10 detector taped to sample in
% epi-direction with a 170µm cover slip and small air-gap.
clear cfg
cfg.nphoton=1e8;
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.issrcfrom0=0;

% Define the material properties
n_518F = 1.5037; % Immersion oil 518F
n_glass = 1.507; % BK7 at 1030nm
n_si = 3.56; % at 1030nm
% 1% Intralipid scattering phantom
%n_scatter = 1.323; %1.323; 1.456
%us_scatter = 2.35; 
%g_scatter = 0.6525; 

% Brain
n_scatter =  1.36;
us_scatter = 50;
g_scatter = 0.9;

cfg.prop  =  [0, 0, 1, 1;            % Zero voxels, air   
              0.001, 0, 1, n_518F;   % Immersion oil 518F
              0.001, 0, 1, n_si;     % Detector sillicone reflection dummy 
              0.001, 0, 1, n_glass;  % 170µm glass cover slip
              0.001, us_scatter, g_scatter, n_scatter]; % Scattering medium
focalPoint_z = 13; % Sample starts at 10, 40µm per step

% define the domain
dimxy = 375; % 15mm in x,y plane
dimz = 250;
cfg.unitinmm = 0.04; % 40µm per voxel
cfg.vol=ones(dimxy,dimxy,dimz); % 15x15x10mm simulation volume
cfg.vol(:,:,1:2) = 2; % Detector silicone
cfg.vol(:,:,3:4) = 1; % Detector / glas gap filled with 518F oil
cfg.vol(:,:,5:9) = 3; % 170µm cover slip
cfg.vol(:,:,10:250) = 4; % Brain

% Cut the 2mm hole in the silicone layer 
r_hole = 1 / cfg.unitinmm; % 1mm detector hole radius
[xi, yi, zi] = meshgrid(1:dimxy, 1:dimxy, 1:2);
dist = sqrt((xi - dimxy/2).^2 + (yi - dimxy/2).^2);
cfg.vol(dist < r_hole) = 1; % Set to air, tag 1

% Isotropic source at focal point to allow for all angles in principle. 
cfg.srctype='isotropic';
cfg.srcpos=[dimxy/2,dimxy/2,focalPoint_z];
cfg.srcdir=[0,0,1];

cfg.savedetflag='dpxv'; % save d: det_id, p: ppath, x: exit position, v: exit direction
cfg.autopilot = 1;

% First, simulate photons reaching the detector from the objective
% Detect photons in epi direction leaving the simulation. Select according
% to real detector dimentions in post selection
% Set boundary conditions per-face (-x,-y,-z,+x,+y,+z) to absorption (a)
% and set face-z as detector (001000).
cfg.bc = 'aaaaaa001000';
cfg.isreflect = 1; % Calculate Frenel reflection
cfg.isspecular = 0;
cfg.maxdetphoton = 1e8;

%mcxpreview(cfg);

% Start the simulation
[flux,detp,vol,seed] = mcxlab(cfg);
 
% Select photons according to dimensions of real detector quadrant for replay
px = (detp.p(:,1)-dimxy/2)*cfg.unitinmm;
py = (detp.p(:,2)-dimxy/2)*cfg.unitinmm;
pz = detp.p(:,3)*cfg.unitinmm;

pd_slit = 0.15; % Tiling slit of detector S11142-10
pd_hole = 2.3;  % Diameter of the hole in detector S11142-10
pd_aa = 14; % Active area size in x and y of S11142-10
% Select photons indecies according to the size of the detector
s1 = px < pd_aa/2 & px > -pd_aa/2 & py < pd_aa/2 & py > -pd_aa/2; % Select photons in square active area
s2 = sqrt(px.^2 + py.^2) > pd_hole/2; % Select photons outside of hole
s3 = (px > pd_slit/2 | px < -pd_slit/2) & (py > pd_slit/2 | py < -pd_slit/2); % Remove photons from tiling slits
sQ12 = py > 0;% & py > 0; % Select the upper two quadrants
% Create a selection mask for left and right side via bitwise AND operation
mask = s1 & s2 & s3 & sQ12;  
clear('s1','s2','s3','sQ12','pd_slit','pd_hole','pd_aa');
% Select photos via mask
px_sel = px(mask);
py_sel = py(mask);
pz_sel = pz(mask);


% Replay with photons reaching the actual detector geometry!
% Modify to detect the replayed photons leaving the focal volume to get
% their directions
cfg_replay = cfg;
cfg_replay.seed = seed.data(:,mask);         % define cfg.seed using the returned seeds
cfg_replay.detphotons= detp.data(:,mask);    % define cfg.detphotons. using the returned detp data
%cfg_replay.outputtype='jacobian';   % tell mcx to output absorption (mu_a) Jacobian

% Sperical zero-shell around illumination point. Set everything 0 except
% for a small sphere
r_0shell = 4; % Radius of zero shell for detection to catch immediate direction of outgoing photons
[xi, yi, zi] = meshgrid(1:dimxy, 1:dimxy, 1:dimz);
dist = sqrt((xi - dimxy/2).^2 + (yi - dimxy/2).^2 + (zi - focalPoint_z).^2);
cfg_replay.vol(dist > r_0shell) = 0; % zero outside of small sphere

% Introduce a new detector to detect the photon directions emitted by the
% source that fell within detector geometry in the previous simulation run.
cfg_replay.detpos=[dimxy/2,dimxy/2,focalPoint_z, r_0shell*2]; % New detector
cfg_replay.bc = 'aaaaaa000000'; % Disable old detector
cfg_replay.isreflect = 0; % Disable reflections
% Also turn of scattering and refraction by setting all µs=0, g=1, and n=1
cfg_replay.prop  =  [0, 0, 1, 1;      % Zero voxels, air   
                     0.001, 0, 1, 1;  % Air
                     0.001, 0, 1, 1;  % Detector sillicone reflection dummy 
                     0.001, 0, 1, 1;  % 170µm glass cover slip
                     0.001, 0, 1, 1]; % Scattring medium

% run replay
[flux_sphRep, detp_sphRep, vol_sphRep, seed_sphRep]=mcxlab(cfg_replay);

% Seperate photons with z+ and z- directions
vz = detp_sphRep.v(:,3); 
s_up = vz > 0; % Select photons going up
s_down = vz < 0; % Select photons going not up

% Copy original setup for replay of where upwards launched photons arrive
cfg_replay_up = cfg; 
cfg_replay_up.seed = seed_sphRep.data(:,s_up);        
cfg_replay_up.detphotons= detp_sphRep.data(:,s_up);

% Copy original setup for replay of where downwards launched photons arrive
cfg_replay_down = cfg;
cfg_replay_down.seed = seed_sphRep.data(:,s_down);        
cfg_replay_down.detphotons= detp_sphRep.data(:,s_down);

[flux_up, detp_up, vol_up, seed_up]=mcxlab(cfg_replay_up);
[flux_down, detp_down, vol_down, seed_down]=mcxlab(cfg_replay_down);

imgSize = round(dimxy/3);
imgMid = round(dimxy/2);
%mcxplotvol(log(flux_up.data(imgMid-imgSize:imgMid+imgSize,imgMid-imgSize:imgMid+imgSize,1:50)));
mcxplotvol(log(flux_up.data));
lighting none
title("Up photons");
%mcxplotvol(log(flux_down.data(imgMid-imgSize:imgMid+imgSize,imgMid-imgSize:imgMid+imgSize,1:50)));
mcxplotvol(log(flux_down.data));
lighting none
title("Down photons");  

fprintf("Photon going up: %d\nPhotons going down: %d\n",numel(detp_up.detid),numel(detp_down.detid))
%% Quick view at the DPLD
v = detp_sphRep.v(detp_sphRep.v(:,3) > 0,:); % Select photons with positive Z direction
MCXData.PhotonData.v = v;

nb = 256;

kBins = linspace(-1,1-2/(nb),nb);

% Corrospondence arrays in units 1 = 2/lambda (here rescaled to one)
uxQ = discretize(v(:,1), kBins);
uyQ = discretize(v(:,2), kBins);

S = zeros(numel(kBins),numel(kBins)); % Source distribution
% Calculate direction distributions from corrospondence arrays
for i = 1:numel(uxQ)
    if(isnan(uxQ(i)) || isnan(uyQ(i)))
    else
        % Obliquity factor for sperical projection in inverse units of
        % 2/lambda (here rescaled to one)
        obFac = sqrt((1)^2 - kBins(uxQ(i))^2 - kBins(uyQ(i))^2);
        obFac = obFac*~imag(obFac); % remove imaginary entries caused by discreet binning
        S(uxQ(i),uyQ(i)) = S(uxQ(i),uyQ(i)) + obFac; % Increase count by +1*obFac
    end
end
S = S + flipud(S)-fliplr(S+flipud(S)); % Construct differential S
S = S/max(S,[],'all');

img = S;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
fprintf("Total contributed Photons: %d \n",size(v,1));
daspect([1 1 1])
colormap jet
colorbar

figure;
plot(kBins,S(128,:));
figform;

%% Save the illumination direction data v as JSON
v = detp_sphRep.v(detp_sphRep.v(:,3) > 0,:); % Select photons with positive Z direction
MCXData.PhotonData.v = v;
savejson('MCXData',MCXData,'FileName','Brain.jdat','Compression','gzip')

%%
figure;
number = 1e6;
px_plot = (detp.p(1e7:1e7+number,1)-dimxy/2)*cfg.unitinmm;
py_plot = (detp.p(1e7:1e7+number,2)-dimxy/2)*cfg.unitinmm;
px_sel_plot = px_plot(mask(1e7:1e7+number));
py_sel_plot = py_plot(mask(1e7:1e7+number));


scatter(px_plot,py_plot,'.');
hold on;
scatter(px_sel_plot,py_sel_plot,'.','r');
xlim([-10,10]);
ylim([-10,10]);
axis equal
%figform;
save_as(gcf,'Scatterplot_detectedSelectedPhotons');

%%

mcxplotvol(log(flux_up.data));
