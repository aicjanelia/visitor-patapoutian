clc
clear
close all

% sourceData will be used to draw reasonable distributions, so pick a
% representatie data set.
sourceData = 'Y:\Rachel\Patapoutian_iPALM_data_analyzed\iPALM_data_processed_EM\CombinedParticles\negYODA1\negYODA1_particles.mat';
% all of the neg YODA1 particles combined
load(sourceData)

%% Get reasonable parameters from the experimental data

%%%%%%%%% BLADE SPACING
closedDist = 20; % nm between closed blades for setting the three locations

% perfect triangle in the xy plane:
xTri = [0,0,0 ; closedDist,0,0 ; 10,sqrt(3/4)*closedDist,0];
xTri(:,1) = xTri(:,1) - mean(xTri(:,1)); % center the triangle for easier rotations
xTri(:,2) = xTri(:,2) - mean(xTri(:,2)); % center the triangle for easier rotations

figure(1)
plot(xTri(:,1),xTri(:,2),'-o')

%%%%%%%%% LOCALIZATION SIGMA
sigma = NaN*ones(8500,2);
c = 0;
for kk = 1:length(subParticles)
    a = subParticles{kk}.sigma;
    b = length(a);
    sigma(c+1:b+c,:) = a;
    c = c+b;
end
% clear kk particles subParticles

%%%% xy sigma
pdxy = fitdist(sigma(:,1),'Kernel');
x = min(sigma(:,1)):.01:max(sigma(:,1));
y = pdf(pdxy,x);

% figure(1)
% histogram(sigma(:,1),'Normalization','pdf')
% hold on
% plot(x,y,'LineWidth',2)
% hold off
% xlabel('\sigma_{x,y} (nm)','FontSize',16)
% set(gca,'FontSize',16)

%%%% z sigma
pdz = fitdist(sigma(:,2),'Kernel');
x = min(sigma(:,2)):.01:max(sigma(:,2));
y = pdf(pdz,x);

% figure(2)
% histogram(sigma(:,2),'Normalization','pdf')
% hold on
% plot(x,y,'LineWidth',2)
% hold off
% xlabel('\sigma_{z} (nm)','FontSize',16)
% set(gca,'FontSize',16)

%%%%%%%%% ROTATION DISTRUBITIONS

pdthetaZ = makedist('Uniform',-pi,pi);
pdthetaXY = makedist('Normal',0,pi/20);

% x = -pi:.01:pi;
% y = pdf(pdthetaXY,x);
% plot(x,y)

% xNow = xTri';
% angX = random(pdthetaXY);
% angY = random(pdthetaXY);
% angZ = random(pdthetaZ);
% 
% xNow = basicRotation(angZ,'Z')*xNow;
% xNow = basicRotation(angY,'Y')*xNow;
% xNow = basicRotation(angX,'X')*xNow;
% 
% xNow = xNow';
% 
% plot3(xTri(:,1),xTri(:,2),xTri(:,3),'-o','LineWidth',3)
% hold on
% plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
% hold off

%%%%%%%%% LOCALIZATION DISTRUBITIONS
L = NaN*ones(length(subParticles),1);
for kk = 1:length(subParticles)
    L(kk) = size(subParticles{kk}.points,1);
end
pdL = fitdist(L(:),'Kernel');
x = min(L(:)):max(L(:));
y = pdf(pdL,x);

% histogram(L,'Normalization','pdf')
% hold on
% plot(x,y,'LineWidth',2)
% hold off

%% Make idealistic PIEZOs
clear particles subParticles

locPerBlade = 5; % localizations per blade
N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['synthetic_' num2str(N) 'particles_' num2str(locPerBlade) 'locPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    xNow = xTri';
    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTri(:,1),xTri(:,2),xTri(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    xNow = repmat(xNow,[locPerBlade 1 1]);
    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make idealistic PIEZOs but asymmetric
clear particles subParticles

% asymmetric triangle in the xy plane:
xTriA = [0,0,0 ; closedDist,0,0 ; 10,sqrt(2)*closedDist,0];
xTriA(:,1) = xTriA(:,1) - mean(xTriA(:,1)); % center the triangle for easier rotations
xTriA(:,2) = xTriA(:,2) - mean(xTriA(:,2)); % center the triangle for easier rotations

locPerBlade = 5; % localizations per blade
N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['syntheticAsymmetric_' num2str(N) 'particles_' num2str(locPerBlade) 'locPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    xNow = xTriA';
    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTriA(:,1),xTriA(:,2),xTriA(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    xNow = repmat(xNow,[locPerBlade 1 1]);
    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make uneven localizations PIEZOs
clear particles subParticles

% locPerBlade = makedist('Normal',5,2); % localizations per blade
locPerBlade = makedist('Normal',5,3); % localizations per blade
N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

% saveTag = ['synthetic_' num2str(N) 'particles_norm5Dist2locPerBlade_scatterSig1and1-5'];
saveTag = ['synthetic_' num2str(N) 'particles_norm5Dist3locPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    xNow = xTri';
    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTri(:,1),xTri(:,2),xTri(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    for jj = 1:size(xNow,1)
        n = round(random(locPerBlade));
        if n > 1 % if 1, no need to repeat, if <1, round up to 1
            xNow = [xNow ; repmat(xNow(jj,:),[n-1 1])];
        end
    end

    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make uneven localization PIEZOs but asymmetric
clear particles subParticles

% asymmetric triangle in the xy plane:
xTriA = [0,0,0 ; closedDist,0,0 ; 10,sqrt(2)*closedDist,0];
xTriA(:,1) = xTriA(:,1) - mean(xTriA(:,1)); % center the triangle for easier rotations
xTriA(:,2) = xTriA(:,2) - mean(xTriA(:,2)); % center the triangle for easier rotations

locPerBlade = makedist('Normal',5,2); % localizations per blade
% locPerBlade = makedist('Normal',5,3); % localizations per blade
N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['syntheticAsymmetric_' num2str(N) 'particles_norm5Dist2locPerBlade_scatterSig1and1-5'];
% saveTag = ['syntheticAsymmetric_' num2str(N) 'particles_norm5Dist3locPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    xNow = xTriA';
    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTriA(:,1),xTriA(:,2),xTriA(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    for jj = 1:size(xNow,1)
        n = round(random(locPerBlade));
        if n > 1 % if 1, no need to repeat, if <1, round up to 1
            xNow = [xNow ; repmat(xNow(jj,:),[n-1 1])];
        end
    end

    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make mixed population of sym & asym, but clear number of localizations per blade
clear particles subParticles

% asymmetric triangle in the xy plane:
xTriA = [0,0,0 ; closedDist,0,0 ; 10,sqrt(2)*closedDist,0];
xTriA(:,1) = xTriA(:,1) - mean(xTriA(:,1)); % center the triangle for easier rotations
xTriA(:,2) = xTriA(:,2) - mean(xTriA(:,2)); % center the triangle for easier rotations

locPerBlade = 5; % localizations per blade
N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['syntheticMixedPopulation_' num2str(N) 'particles_' num2str(locPerBlade) 'locPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    if mod(kk,2)
        xNow = xTri';
    else
        xNow = xTriA';
    end

    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTriA(:,1),xTriA(:,2),xTriA(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    xNow = repmat(xNow,[locPerBlade 1 1]);
    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make mixed population of sym & asym
clear particles subParticles

% asymmetric triangle in the xy plane:
xTriA = [0,0,0 ; closedDist,0,0 ; 10,sqrt(2)*closedDist,0];
xTriA(:,1) = xTriA(:,1) - mean(xTriA(:,1)); % center the triangle for easier rotations
xTriA(:,2) = xTriA(:,2) - mean(xTriA(:,2)); % center the triangle for easier rotations

% locPerBlade = makedist('Normal',5,2); % localizations per blade
locPerBlade = makedist('Normal',5,3); % localizations per blade
N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

% saveTag = ['syntheticMixedPopulation_' num2str(N) 'particles_norm5Dist2locPerBlade_scatterSig1and1-5'];
saveTag = ['syntheticMixedPopulation_' num2str(N) 'particles_norm5Dist3locPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    if mod(kk,2)
        xNow = xTri';
    else
        xNow = xTriA';
    end

    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTriA(:,1),xTriA(:,2),xTriA(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    for jj = 1:size(xNow,1)
        n = round(random(locPerBlade));
        if n > 1 % if 1, no need to repeat, if <1, round up to 1
            xNow = [xNow ; repmat(xNow(jj,:),[n-1 1])];
        end
    end

    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make idealistic barbell (test of the rotation properties)
clear particles subParticles

xDi = [-closedDist/2,0,0 ; closedDist/2,0,0];

locPerBlade = 5; % localizations per blade
N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['syntheticBarbell_' num2str(N) 'particles_' num2str(locPerBlade) 'locPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    xNow = xDi';
    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xDi(:,1),xDi(:,2),xDi(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    xNow = repmat(xNow,[locPerBlade 1 1]);
    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make uneven localizations PIEZOs - draw from data distribution
clear particles subParticles

N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['synthetic_' num2str(N) 'particles_dataDistThirdlocPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    xNow = xTri';
    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTri(:,1),xTri(:,2),xTri(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    for jj = 1:size(xNow,1)
        n = round(random(pdL)/3);
        if n > 1 % if 1, no need to repeat, if <1, round up to 1
            xNow = [xNow ; repmat(xNow(jj,:),[n-1 1])];
        end
    end

    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make uneven localization PIEZOs but asymmetric  - draw from data distribution
clear particles subParticles

% asymmetric triangle in the xy plane:
xTriA = [0,0,0 ; closedDist,0,0 ; 10,sqrt(2)*closedDist,0];
xTriA(:,1) = xTriA(:,1) - mean(xTriA(:,1)); % center the triangle for easier rotations
xTriA(:,2) = xTriA(:,2) - mean(xTriA(:,2)); % center the triangle for easier rotations

N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['syntheticAsymmetric_' num2str(N) 'particles_dataDistThirdlocPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    xNow = xTriA';
    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTriA(:,1),xTriA(:,2),xTriA(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    for jj = 1:size(xNow,1)
        n = round(random(pdL)/3);
        if n > 1 % if 1, no need to repeat, if <1, round up to 1
            xNow = [xNow ; repmat(xNow(jj,:),[n-1 1])];
        end
    end

    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make mixed population of sym & asym  - draw from data distribution
clear particles subParticles

% asymmetric triangle in the xy plane:
xTriA = [0,0,0 ; closedDist,0,0 ; 10,sqrt(2)*closedDist,0];
xTriA(:,1) = xTriA(:,1) - mean(xTriA(:,1)); % center the triangle for easier rotations
xTriA(:,2) = xTriA(:,2) - mean(xTriA(:,2)); % center the triangle for easier rotations

N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['syntheticMixedPopulation_' num2str(N) 'particles_dataDistThirdlocPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    if mod(kk,2)
        xNow = xTri';
    else
        xNow = xTriA';
    end

    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTriA(:,1),xTriA(:,2),xTriA(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    for jj = 1:size(xNow,1)
        n = round(random(pdL)/3);
        if n > 1 % if 1, no need to repeat, if <1, round up to 1
            xNow = [xNow ; repmat(xNow(jj,:),[n-1 1])];
        end
    end

    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')

%% Make mixed population of sym sizes  - draw locPerBlade from data distribution
clear particles subParticles

closedDist2 = 30; % nm between closed blades for setting the three locations

% perfect triangle in the xy plane:
xTri2 = [0,0,0 ; closedDist2,0,0 ; 10,sqrt(3/4)*closedDist2,0];
xTri2(:,1) = xTri2(:,1) - mean(xTri2(:,1)); % center the triangle for easier rotations
xTri2(:,2) = xTri2(:,2) - mean(xTri2(:,2)); % center the triangle for easier rotations

N = 100; % how many particles to make

%%% Make the localization scatter pretty limited for this 'idealized' case
pdxyLoc = makedist('Normal',0,1);
pdzLoc = makedist('Normal',0,1.5);

saveTag = ['syntheticMixedSize_' num2str(N) 'particles_dataDistThirdlocPerBlade_scatterSig1and1-5'];
if ~exist(saveTag,'file')
    mkdir(saveTag)
end

figure(3)
set(gcf,'Position',[550 150 1500 1000])

%%%%% Create the particles
particles = cell(1,N);
subParticles = particles;
for kk = 1:N

    if mod(kk,2)
        xNow = xTri';
    else
        xNow = xTri2';
    end

    angX = random(pdthetaXY);
    angY = random(pdthetaXY);
    angZ = random(pdthetaZ);
    
    xNow = basicRotation(angZ,'Z')*xNow;
    xNow = basicRotation(angY,'Y')*xNow;
    xNow = basicRotation(angX,'X')*xNow;
    
    xNow = xNow';
    
%     plot3(xTriA(:,1),xTriA(:,2),xTriA(:,3),'-o','LineWidth',1)
%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'-s','LineWidth',3)
%     hold off

    for jj = 1:size(xNow,1)
        n = round(random(pdL)/3);
        if n > 1 % if 1, no need to repeat, if <1, round up to 1
            xNow = [xNow ; repmat(xNow(jj,:),[n-1 1])];
        end
    end

    sig = [random(pdxy,[size(xNow,1) 1]) random(pdz,[size(xNow,1) 1])];
    jitt = [random(pdxyLoc,[size(xNow,1) 2]) random(pdzLoc,[size(xNow,1) 1])];
    xNow = xNow+jitt;

%     hold on
%     plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
%     hold off

    % Put into the Heydarian "subParticle" format
    subParticles{1,kk}.points = xNow;
    subParticles{1,kk}.sigma = sig;
    subParticles{1,kk}.txtData = NaN;
    subParticles{1,kk}.headers = NaN;

    % Put into the Heydarian "particles" format
    particles{1,kk}.coords(:,1:3) = xNow;
    particles{1,kk}.coords(:,5) = sig(:,1);
    particles{1,kk}.coords(:,10) = sig(:,2);

    if kk <= 25
        subplot(5,5,kk)
        plot3(xNow(:,1),xNow(:,2),xNow(:,3),'*','LineWidth',3)
        set(gca,'DataAspectRatio',[1 1 1])
    end

end

save([saveTag filesep saveTag '_particles.mat'],'particles','subParticles',...
    'pdxyLoc','pdzLoc','pdthetaXY','pdthetaZ','pdxy','pdz')
saveas(gcf,[saveTag filesep saveTag '_exampleParticles.png'],'png')


%%
