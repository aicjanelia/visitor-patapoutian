%% INSTRUCTIONS
% (1) Run beadRemoval_v0_anisotropic to remove beads from txt data
% (2) Load bead-removed txt data back into PeakSelector
% (3) Select option: PeakSelector > SpecialFunctions > Swap Z with Unwrapped Z
% (4) Open the "Cust. Tiff" window.
%       Make sure "Filter" is set to "Frame Peaks"
%       Set "NM per Image Pixel" to 2
% (5) Click "Render" --> this step could be slow
% (6) Click "Save TIFF float" and save file
% (7) Fill in USER PARAMETERS in the script below, run piezoSegment_beadsRemoved
%       Note: want the _r file from rendering (red channel)

clc, clear, close all
%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Filenames
dataDir = 'Z:\Rachel\Patapoutian\21.11.02-5\Run1-561\beadRemoval_piezoSegment_averaged\'; % Where is the data?
saveTag = 'Run1-561_c123_sum_X14_processed_overlay_Fiducial_transform_complete_IDL_ASCII_200-400unwZ_beadsRemoved_rRemoveX3_rRemoveY7'; % Base filename for the data
% dataDir & saveTag will be used in the file name used when saving figures, etc.

%%%%% Localization Data
filenameTxt = [dataDir saveTag '.txt']; % Filename (including path if not on current MATLAB path) for the bead-removed txt data
xyScale = 133.33; % nm/pix in the txt file (which is different than the rendered image!)

%%%%% Rendered Image
filenameImg = [dataDir saveTag '_r.tif']; % Filename (including path if not on current MATLAB path) for the rendered image
nmPix = 3; % nm per pixel in the rendered image

%%%%% Rendered Peak Finding
rBlade = 20/nmPix; % pixels, approximate size of a blade-blob
peakThresh = 0.001; % threshold on the bandpassed image
rThresh = rBlade/2; % pixels, don't allow two peaks within the size of the blade

%%%%% Neighbor requirements
neighborNumber = 2; % Number of neighbors each peak should have
neighborDist = 60; % nm, maximum neighbor distance
minDist = 9; % nm; zero gets rid of self comparisons, can also set higher to get rid of artifacts
bord = 40/nmPix; % pixels, border to add around the candidates during segmentation

%% First find peaks in the rendered data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imR = imread(filenameImg,'tif');
imR = double(imR);

% imagesc(R)
% set(gca,'DataAspectRatio',[1 1 1])
% caxis([0 .1])

%%%% bpass rendered data
bR = bpass(imR,0,rBlade); % 0 = no smoothing

% imagesc(bR)
% set(gca,'DataAspectRatio',[1 1 1])
% caxis([0 .1])

%%%% pknd rendered data
pkR = pkfnd(bR,peakThresh,rThresh);

% Check on the peak finding
figure(101)
imagesc(imR)
set(gca,'DataAspectRatio',[1 1 1])
caxis([0 .1])
hold on
plot(pkR(:,1),pkR(:,2),'om','MarkerFaceColor','m')
hold off
title('Sanity Check: Are peaks being found and not noise?')

%% Set requirements on neighboring peaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = squareform(pdist(pkR)); % All distances between neighbors
Zplus2 = (Z < neighborDist/nmPix) & (Z > minDist/nmPix);

% Set a requirment on neighbor number
% plus2 = sum(Zplus2,2) >= 2; % This would keep cases where 2 piezos are close together
plus2 = sum(Zplus2,2) == neighborNumber; % Stricter requirement = cleaner localization segmentation???

figure(102)
imagesc(imR)
set(gca,'DataAspectRatio',[1 1 1])
caxis([0 .1])
hold on
plot(pkR(:,1),pkR(:,2),'om','MarkerFaceColor','m','MarkerSize',2)
plot(pkR(plus2,1),pkR(plus2,2),'sg','LineWidth',2)
hold off
title('Sanity Check: Does the neighborfinding remove too many or not enough points?')

%% Segement points with correct number of neighbors

ids = find(plus2);
count = 1;
lims = NaN*zeros(1,4);

% loop through those points that meet plus2 criteria
while ~isempty(ids)

    % rows and columns of Z are the same
    idNow = ids(1);
    neigh = find(Zplus2(idNow,:));
    if length(neigh) ~= neighborNumber
        error('mismatch in neighbor number')
    end

    % Stricter condition:
    % all must be part of plus2 AND _only_ neighbor each other
    n1 = sort([idNow neigh]);
    n2 = sort([neigh(1) find(Zplus2(neigh(1),:))]);
    n3 = sort([neigh(2) find(Zplus2(neigh(2),:))]);

    if isequal(n1,n2) && isequal(n1,n3)
        x = pkR([idNow neigh],1);
        y = pkR([idNow neigh],2);
        lims(count,:) = [min(x(:))-bord max(x(:))+bord min(y(:))-bord max(y(:))+bord];
        count = count+1;
    end

    % remove the ids that were used from the list
    ids(1) = []; % current particle
    for kk = 1:neighborNumber % already used neighbors
        ids(ids == neigh(kk)) = [];
    end

end
clear n1 n2 n3 neigh count ids idNow x y

%% Rendered images of first 25 of the segmented points
Nfig = min(length(lims),25);
figure(103)
for kk = 1:Nfig

    subplot(5,5,kk)
    Rslice = imR(:,lims(kk,1):lims(kk,2));
    Rslice = Rslice(lims(kk,3):lims(kk,4),:);
    imagesc(Rslice)
    set(gca,'DataAspectRatio',[1 1 1])
    caxis([0 .25])
    caxis([0 0.04])

end
clear Rslice

%% Connect Segementation to Localizations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load ASCII data
txtData = readmatrix(filenameTxt);
headers = detectImportOptions(filenameTxt);
headers = headers.VariableNames;
% Format the data appropriately
x = txtData(:,3)*xyScale;
y = txtData(:,4)*xyScale;
y = (size(imR,1)+1)*nmPix - y; % flip the coordinate system to match the image...
z = txtData(:,45);
sigXY = xyScale*(txtData(:,15)+txtData(:,16))/2;
sigZ = txtData(:,36);
blades = txtData(:,27) == 1; % blades are the red channel, may look at central pores later

% Visualize localizations on image (check for rough registration)
figure(104)
imagesc(imR)
caxis([0 0.001])
set(gca,'DataAspectRatio',[1 1 1])
hold on
plot(x(blades)/nmPix,y(blades)/nmPix,'.m')
hold off
title('Sanity Check: Do the registrations roughly register to the rendered image?')

%% Find the localizations that match the segmented candidate images

particles = cell(1,length(lims));
subParticles = particles;
lims2 = lims*nmPix;

for kk = 1:length(lims)

    rowsNowAll = x>lims2(kk,1) & x<lims2(kk,2) & y>lims2(kk,3) & y<lims2(kk,4); % spatial location
    rowsNow = rowsNowAll & blades; % red channel <-- might revist other channel later, stay simple for now
    if sum(rowsNow) == 0
        error('No particles') % Unlikely error, but might use this to set other conditions later...
    end

    % Center the particle in its own reference frame
    xnow = x(rowsNow);
    xbar = mean(xnow);
    xnow = xnow-xbar;

    ynow = y(rowsNow);
    ybar = mean(ynow);
    ynow = ynow-ybar;

    znow = z(rowsNow);
    zbar = mean(znow);
    znow = znow-zbar;

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = [xnow ynow znow];
    subParticles{1,kk}.sigma = [sigXY(rowsNow) sigZ(rowsNow)];
    % Have all the peakselector data tag along
    subParticles{1,kk}.txtData = txtData(rowsNow,:);
    subParticles{1,kk}.headers = headers;
    
    % This is for comparing to the segmented image coordinate system
    xnow = x(rowsNow)-lims2(kk,1);
    ynow = y(rowsNow)-lims2(kk,3);
    znow = z(rowsNow) - zbar;
    subParticles{1,kk}.image = [xnow ynow znow];

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = subParticles{1,kk}.points;
    particles{1,kk}.coords(:,5) = sigXY(rowsNow);
    particles{1,kk}.coords(:,10) = sigZ(rowsNow);

    % Have the second channel tag along in case we want to use this information
    rowsNow2 = rowsNowAll & ~blades;
    xnow = x(rowsNow2) - xbar;
    ynow = y(rowsNow2) - ybar;
    znow = z(rowsNow2) - zbar;
    subParticles{1,kk}.pointsCh2 = [xnow ynow znow];
    subParticles{1,kk}.sigmaCh2 = [sigXY(rowsNow2) sigZ(rowsNow2)];
    subParticles{1,kk}.txtDataCh2 = txtData(rowsNow2,:);
    % For image comparison with second channel
    xnow = x(rowsNow2)-lims2(kk,1);
    ynow = y(rowsNow2)-lims2(kk,3);
    znow = z(rowsNow2) - zbar;
    subParticles{1,kk}.imageCh2 = [xnow ynow znow];


end
clear rowsNow xnow ynow znow

%% Rendered images of first 25 of the segmented points
% + the localization data
figure(105)
for kk = 1:Nfig
    
    Rslice = imR(:,lims(kk,1):lims(kk,2));
    Rslice = Rslice(lims(kk,3):lims(kk,4),:);

    xyz = subParticles{1,kk}.image/nmPix;
    xyzCh2 = subParticles{1,kk}.imageCh2/nmPix;

    subplot(5,5,kk)
    imagesc(Rslice)
    set(gca,'DataAspectRatio',[1 1 1])
    caxis([0 0.04]) % optimized for the 2 nm rendering
    colormap gray
    hold on
        plot(xyz(:,1),xyz(:,2),'.m','MarkerSize',8)
        plot(xyzCh2(:,1),xyzCh2(:,2),'.g','MarkerSize',8)
    hold off

end 
clear Rslice kk

%% Save Segmentation Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([dataDir saveTag '_segmentationWorkspace.mat'],'-v7.3') % This is a large file with _everything_
save([dataDir saveTag '_particles.mat'],'particles','subParticles') % This is a small file that can be used for downstream particle averaging