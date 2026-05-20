function S = calcSource(lambda,FOVx,pxN,sim_data_path)

% Calculate lateral frequency grids
dx = FOVx/pxN; % Image pixel size in µm
kx_max = 0.5/dx; % Max spatial frequency in lines/µm after fft of image
dkx = 2*kx_max/pxN; % Frequency step per pixel in lines/µm of image fft
kxBins = linspace(-kx_max,kx_max-dkx,pxN); % image frequency space k axis

% Load simulation file to calculate source function S
sample_sim = loadjson(sim_data_path);
fprintf("Simulation data loaded!\nPath: %s\n",sim_data_path);
% Get data out of simulation file
v = sample_sim.MCXData.PhotonData.v;
% Corrospondence arrays in units 1 = 2/lambda
uxQ = discretize(2/lambda*v(:,1), kxBins);
uyQ = discretize(2/lambda*v(:,2), kxBins);

S = zeros(numel(kxBins),numel(kxBins)); % Source distribution
% Calculate direction distributions from corrospondence arrays
for i = 1:numel(uxQ)
    if(isnan(uxQ(i)) || isnan(uyQ(i)))
    else
        % Obliquity factor for sperical projection in inverse units of lambda
        obFac = sqrt((2/lambda)^2 - kxBins(uxQ(i))^2 - kxBins(uyQ(i))^2);
        obFac = obFac*~imag(obFac); % remove imaginary entries caused by discreet binning
        S(uxQ(i),uyQ(i)) = S(uxQ(i),uyQ(i)) + obFac; % Increase count by +1*obFac
    end
end
S = (S + flipud(S)) - fliplr(S + flipud(S)); % Construct differential S from all four illumination directions
S = S/max(S,[],'all');

end