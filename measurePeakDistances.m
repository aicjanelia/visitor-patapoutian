%%%%%% INSTRUCTIONS
% (1) Run batch_piezoSymFold over the data you wish to analyze. The
%       relevant file created by that process is
%       '*_piezoAveragingWorkspaceRotated.mat'
% (2) Fill in the user parameters below and run the script. The variable
%       directories allows this script to loop through multiple experiments,
%       analyzing each one in turn.

clc, clear, close all

%% USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directories = {'Y:\Rachel\Patapoutian_iPALM_data_analyzed\iPALM_data_processed_EM\CombinedParticles\negYODA1\',...
    'Y:\Rachel\Patapoutian_iPALM_data_analyzed\iPALM_data_processed_EM\CombinedParticles\posYODA1_2\'};

scale = 5; % This is the scale used in the averaging process

%% MAKE MEASUREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for mm = 1:length(directories)

    directory = directories{mm};
    fileName = dir([directory '*_piezoAveragingWorkspaceRotated.mat']);
    if isempty(fileName)
        error('Rotated SuperParticle Not Found')
    end
    fileName = fileName(1).name;
    saveTag = fileName(1:end-35);

    %% Load the raw data
    
    load([directory fileName])

    data = superParticleWithPK3Rot{1,6};    
    density = mvksdensity(data,data,'Bandwidth',scale/2); %  This is the probability density
    
    %% Calculate a threshold
    
    h = histogram(density(:));
    T = otsuthresh(h.BinCounts);
    blades = density > T*max(density(:));

    x = data(blades,1);
    y = data(blades,2);
    z = data(blades,3);
    d = density(blades);
    
    figure(1)
    set(gcf,'Position',[500 275 560*2 420*2])    
    scatter3(x,y,z, d*5*10^6, d, '.')
    colormap(cool)
    xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'FontSize',16)
    colorbar       
    set(gcf,'Color','white')
    
    view([46 49])
    xlim(25*[-1 1])
    ylim(25*[-1 1])
    zlim(10*[-1 1])
    caxis((10^-5)*[2 6])
    
    saveas(gcf,[directory saveTag '_OtsuSuperParticle.png'],'png')
    saveas(gcf,[directory saveTag '_OtsuSuperParticle.fig'],'fig')
    
    %% Create Gridded Data to Calculate Coarse Peaks
    
    meshPos = ceil(max([max(x) abs(min(x)) max(y) abs(min(y))]))+5; % Max position on the grid
    meshSpace = 1; % 1 is so indices = nm (could consider other options if necessary)
    
    [Xq,Yq,Zq] = meshgrid(-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos);
    Vq = mvksdensity(data,[Xq(:) Yq(:) Zq(:)],'Bandwidth',scale/2);    
    densityGrid = reshape(Vq,size(Xq));

    figure(2)
    imagesc(densityGrid(:,:,meshPos))
    set(gca,'DataAspectRatio',[1 1 1])
    saveas(gcf,[directory saveTag '_densityGrid_slice' num2str(meshPos) '.png'],'png')
    saveas(gcf,[directory saveTag '_densityGrid_slice' num2str(meshPos) '.fig'],'fig')
    
%     h = histogram(densityGrid(:));
%     T2 = otsuthresh(h.BinCounts);    
%     T2 = T2*max(densityGrid(:));
    T2 = T*max(density(:));

    % Find pixels that are above the threshold
    k = find(densityGrid>T2);
    [a,b,c]=ind2sub(size(densityGrid),k);
    edgeSize = max(size(densityGrid));

    keepers = 0*a; % Values that count as peaks
    for kk = 1:length(a)
        x1=a(kk);
        x2=b(kk);
        x3=c(kk);
        % This is an adaptation of the algorithm in pkfnd, adapted to 3D data
        if x1>1 && x1<edgeSize && x2>1 && x2<edgeSize && x3>1 && x3<edgeSize % avoid the image edges        
            if densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2,x3) && ... % all 26 neighbors
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2-1,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2+1,x3) && ...    
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2-1,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2+1,x3) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2-1,x3) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2+1,x3) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2-1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2-1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1,x2+1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1,x2+1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2-1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2-1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2+1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1-1,x2+1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2-1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2-1,x3+1) && ...
               densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2+1,x3-1) && densityGrid(x1,x2,x3) >= densityGrid(x1+1,x2+1,x3+1)
    
                    keepers(kk) = 1; % If the pixel is the brightest of all its neighbors, keep it
    
            end
    
        end
    end
    
    pkx = [a(keepers>0) b(keepers>0) c(keepers>0)]; % Peaks in the mesh grid coordinate system
    
    %% Calculate Fine Peaks Based on Point Cloud 

    sz = scale; %window for fitting neighbor peaks

    pkC = NaN*ones(size(pkx,1),4);
    for kk = 1:sum(keepers)

        r = sqrt(sum((data(:,1)-pkx(kk,2)+meshPos).^2 + (data(:,2)-pkx(kk,1)+meshPos).^2 + (data(:,3)-pkx(kk,3)+meshPos).^2,2));
        neigh = (r<sz);
        if ~sum(neigh)
            error('No neighbors found around peak')
        end

        tmp = density(neigh);
        norm = sum(tmp);
        xavg = sum(data(neigh,1).*tmp)./norm; % Weighted average based on local density
        yavg = sum(data(neigh,2).*tmp)./norm;
        zavg = sum(data(neigh,3).*tmp)./norm;

        %  Radius of gyration (= sum(m_i * r_i^2)/sum(m_i)), m = mass, r = radius
        r = sqrt(sum((data(neigh,1)-xavg).^2 + (data(neigh,2)-yavg).^2 + (data(neigh,3)-zavg).^2,2));
        rg = sum(tmp.*r.^2)/norm;

        pkC(kk,:) = [xavg yavg zavg rg]; % This is in nm relative in the superparticle coordinate system

    end
%     disp(pkC)

    %% Plot locations of found peaks
    
    toPlotCoarse = [b(keepers>0) a(keepers>0) c(keepers>0)]-meshPos; % Move back to particle coordinate system, row = y, col = x
    toPlotFine = [pkC(:,1) pkC(:,2) pkC(:,3)];
    
    figure(3)
    set(gcf,'Position',[500 275 560*2 420*2])    
    scatter3(x,y,z, d*5*10^6, d, '.')
    hold on
    l1 = plot3(toPlotCoarse(:,1),toPlotCoarse(:,2),toPlotCoarse(:,3),'ok','MarkerSize',10,'LineWidth',2);
    l2 = plot3(toPlotFine(:,1),toPlotFine(:,2),toPlotFine(:,3),'k^','MarkerSize',10,'LineWidth',2);
    hold off
    colormap(cool)
    xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
    set(gca,'DataAspectRatio',[1 1 1])
    set(gca,'FontSize',16)
    colorbar       
    set(gcf,'Color','white')
    legend([l1 l2],'Coarse Peaks','Fine Peaks')
    
    view(2)
    
    saveas(gcf,[directory saveTag '_FoundPeaksOnParticle.png'],'png')
    saveas(gcf,[directory saveTag '_FoundPeaksOnParticle.fig'],'fig')

    dPeak = mvksdensity(data,pkC(:,1:3),'Bandwidth',scale/2); %  This is the probability density ath the peaks
    dPeakCoarse = mvksdensity(data,pkC(:,1:3),'Bandwidth',scale/2); %  This is the probability density ath the peaks

    figure(4)
    scatter3(data(:,1),data(:,2),density(:),[],density(:))
    set(gca,'DataAspectRatio',[1 1 10^-6])
    hold on
    for kk = 1:size(pkC,1)
        l1 = plot3(toPlotCoarse(kk,1),toPlotCoarse(kk,2),dPeakCoarse(kk),'ok','MarkerSize',10);
        l2 = plot3(pkC(kk,1),pkC(kk,2),dPeak(kk),'pr','MarkerFaceColor','r','MarkerSize',10);
    end
    hold off
    legend([l1 l2],'Coarse Peaks','Fine Peaks')
    saveas(gcf,[directory saveTag '_FoundPeaks2D.png'],'png')
    saveas(gcf,[directory saveTag '_FoundPeaks2D.fig'],'fig')

    %% Intermediate clean up
    close all

    %% Bootstrap a confidence interval for the peak locations
    disp('Starting Bootstrap')

    s = RandStream('mlfg6331_64'); % Reproducible random stream
    Nresample = 1000;
    indSample = 1:size(data,1);
    pkCBoot = cell(Nresample,1);
    dataBoot = pkCBoot; densityBoot = pkCBoot; bladesBoot = pkCBoot;
    xBoot = pkCBoot; yBoot = xBoot; zBoot = xBoot; dBoot = xBoot;
    dPeakBoot = pkCBoot;

    for ii = 1:Nresample
        disp(ii)
%         tic

        %%% Randomly sample with replacement to a dataset of the same size
        dataSamp = datasample(s,indSample,size(data,1),'Replace',true);
        dataSamp = data(dataSamp,:);
        
        densitySamp = mvksdensity(dataSamp,dataSamp,'Bandwidth',scale/2); %  This is the probability density
    
        %%% Calculate a threshold        
        figure(10)
        h = histogram(densitySamp(:));
        T = otsuthresh(h.BinCounts);
        bladesSamp = densitySamp > T*max(densitySamp(:));
    
        xSamp = dataSamp(bladesSamp,1);
        ySamp = dataSamp(bladesSamp,2);
        zSamp = dataSamp(bladesSamp,3);
        dSamp = densitySamp(bladesSamp);
        
%         figure(1)
%         set(gcf,'Position',[500 275 560*2 420*2])    
%         scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
%         colormap(cool)
%         xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
%         set(gca,'DataAspectRatio',[1 1 1])
%         set(gca,'FontSize',16)
%         colorbar       
%         set(gcf,'Color','white')
%         drawnow
        
        view([46 49])
        xlim(25*[-1 1])
        ylim(25*[-1 1])
        zlim(10*[-1 1])
        caxis((10^-5)*[2 6])
    
        %%% Create Gridded Data to Calculate Coarse Peaks        
        meshPos = ceil(max([max(xSamp) abs(min(xSamp)) max(ySamp) abs(min(ySamp))]))+5; % Max position on the grid
        meshSpace = 1; % 1 is so indices = nm (could consider other options if necessary)
        
        [Xq,Yq,Zq] = meshgrid(-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos);
        Vq = mvksdensity(dataSamp,[Xq(:) Yq(:) Zq(:)],'Bandwidth',scale/2);    
        densityGridSamp = reshape(Vq,size(Xq));
    
%         figure(2)
%         imagesc(densityGridSamp(:,:,meshPos))
%         set(gca,'DataAspectRatio',[1 1 1])
        
        figure(10)
        h = histogram(densityGridSamp(:));
        T2 = otsuthresh(h.BinCounts);    
        T2 = T2*max(densityGridSamp(:));
        T2 = T*max(density(:));

        % Don't allow isolated components, only want to keep the main PIEZO
        % structure, not small peaks out in the 'noise'
        vol = densityGridSamp>T2;
        CC = bwconncomp(vol);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        filtered_vol = false(size(vol));
        filtered_vol(CC.PixelIdxList{idx}) = true;
    
        % Find pixels that are above the threshold
        k = find(filtered_vol);
        [a,b,c]=ind2sub(size(densityGridSamp),k);
        edgeSize = max(size(densityGridSamp));
    
        keepers = 0*a; % Values that count as peaks
        for kk = 1:length(a)
            x1=a(kk);
            x2=b(kk);
            x3=c(kk);
            % This is an adaptation of the algorithm in pkfnd, adapted to 3D data
            if x1>1 && x1<edgeSize && x2>1 && x2<edgeSize && x3>1 && x3<edgeSize % avoid the image edges        
                if densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3) && ... % all 26 neighbors
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3) && ...    
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3+1)
        
                        keepers(kk) = 1; % If the pixel is the brightest of all its neighbors, keep it
        
                end
        
            end
        end
        
        pkxSamp = [a(keepers>0) b(keepers>0) c(keepers>0)]; % Peaks in the mesh grid coordinate system
        
        %%% Calculate Fine Peaks Based on Point Cloud 
    
        sz = scale; %window for fitting neighbor peaks
    
        pkCSamp = NaN*ones(size(pkxSamp,1),4);
        for kk = 1:sum(keepers)
    
            r = sqrt(sum((dataSamp(:,1)-pkxSamp(kk,2)+meshPos).^2 + (dataSamp(:,2)-pkxSamp(kk,1)+meshPos).^2 + (dataSamp(:,3)-pkxSamp(kk,3)+meshPos).^2,2));
            neigh = (r<sz);
            if ~sum(neigh)
                error('No neighbors found around peak')
            end
    
            tmp = densitySamp(neigh);
            norm = sum(tmp);
            xavg = sum(dataSamp(neigh,1).*tmp)./norm; % Weighted average based on local density
            yavg = sum(dataSamp(neigh,2).*tmp)./norm;
            zavg = sum(dataSamp(neigh,3).*tmp)./norm;
    
            %  Radius of gyration (= sum(m_i * r_i^2)/sum(m_i)), m = mass, r = radius
            r = sqrt(sum((dataSamp(neigh,1)-xavg).^2 + (dataSamp(neigh,2)-yavg).^2 + (dataSamp(neigh,3)-zavg).^2,2));
            rg = sum(tmp.*r.^2)/norm;
    
            pkCSamp(kk,:) = [xavg yavg zavg rg]; % This is in nm relative in the superparticle coordinate system
    
        end
%         disp(pkCSamp)

        figure(3)
        set(gcf,'Position',[500 275 560*2 420*2])    
        scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
        hold on
        plot3(pkCSamp(:,1),pkCSamp(:,2),pkCSamp(:,3),'pr','MarkerFaceColor','r','MarkerSize',10);
        hold off
        colormap(cool)
        xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
        set(gca,'DataAspectRatio',[1 1 1])
        set(gca,'FontSize',16)
        colorbar       
        set(gcf,'Color','white')
        view(2)
        drawnow

        dPeakSamp = mvksdensity(dataSamp,pkCSamp(:,1:3),'Bandwidth',scale/2); %  This is the probability density at the peaks
        figure(4)
        scatter3(dataSamp(:,1),dataSamp(:,2),densitySamp(:),[],densitySamp(:))
        set(gca,'DataAspectRatio',[1 1 10^-6])
        hold on
        for kk = 1:size(pkCSamp,1)
            plot3(pkCSamp(kk,1),pkCSamp(kk,2),dPeakSamp(kk),'pr','MarkerFaceColor','r','MarkerSize',10);
        end
        hold off
        drawnow

        pkCBoot{ii} = pkCSamp;

        dataBoot{ii} = dataSamp;
        densityBoot{ii} = densitySamp;
        bladesBoot{ii} = bladesSamp;
        xBoot{ii} = xSamp;
        yBoot{ii} = ySamp;
        zBoot{ii} = zSamp;
        dBoot{ii} = dSamp;
        dPeakBoot{ii} = dPeakSamp;
        
%         toc

    end

    %% Save the data    
    save([directory saveTag '_PeakFinding.mat'],'data','density','blades','x','y','z','d',...
        'densityGrid','pkx','pkC','toPlotCoarse','toPlotFine','dPeak','dPeakCoarse',...
        'pkCBoot','dataBoot','densityBoot','bladesBoot','xBoot','yBoot','zBoot','dBoot','dPeakBoot')

    %% Leave One Out Analysis...
    
    disp('Starting LOO')
tic
    Nloo = size(data,1);
    pkCLoo = cell(Nloo,1);
    dataLoo = pkCLoo; densityLoo = pkCLoo; bladesLoo = pkCLoo;
    xLoo = pkCLoo; yLoo = xLoo; zLoo = xLoo; dLoo = xLoo;
    dPeakLoo = pkCLoo;

    for ii = 2197:Nloo
        if ~mod(ii,50)
            disp(ii)
        end
%         tic

        % Leave one out
        dataSamp = ones(Nloo,1);
        dataSamp(ii) = 0;
        dataSamp = data(dataSamp>0,:);
        
        densitySamp = mvksdensity(dataSamp,dataSamp,'Bandwidth',scale/2); %  This is the probability density
    
        %%% Calculate a threshold        
        figure(10)
        h = histogram(densitySamp(:));
        T = otsuthresh(h.BinCounts);
        bladesSamp = densitySamp > T*max(densitySamp(:));
    
        xSamp = dataSamp(bladesSamp,1);
        ySamp = dataSamp(bladesSamp,2);
        zSamp = dataSamp(bladesSamp,3);
        dSamp = densitySamp(bladesSamp);
        
%         figure(1)
%         set(gcf,'Position',[500 275 560*2 420*2])    
%         scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
%         colormap(cool)
%         xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
%         set(gca,'DataAspectRatio',[1 1 1])
%         set(gca,'FontSize',16)
%         colorbar       
%         set(gcf,'Color','white')
%         drawnow
%         
%         view([46 49])
%         xlim(25*[-1 1])
%         ylim(25*[-1 1])
%         zlim(10*[-1 1])
%         caxis((10^-5)*[2 6])
    
        %%% Create Gridded Data to Calculate Coarse Peaks        
        meshPos = ceil(max([max(xSamp) abs(min(xSamp)) max(ySamp) abs(min(ySamp))]))+5; % Max position on the grid
        meshSpace = 1; % 1 is so indices = nm (could consider other options if necessary)
        
        [Xq,Yq,Zq] = meshgrid(-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos,-meshPos:meshSpace:meshPos);
        Vq = mvksdensity(dataSamp,[Xq(:) Yq(:) Zq(:)],'Bandwidth',scale/2);    
        densityGridSamp = reshape(Vq,size(Xq));
    
%         figure(2)
%         imagesc(densityGridSamp(:,:,meshPos))
%         set(gca,'DataAspectRatio',[1 1 1])
        
%         figure(10)
%         h = histogram(densityGridSamp(:));
%         T2 = otsuthresh(h.BinCounts);    
%         T2 = T2*max(densityGridSamp(:));
        T2 = T*max(density(:));

        % Don't allow isolated components, only want to keep the main PIEZO
        % structure, not small peaks out in the 'noise'
        vol = densityGridSamp>T2;
        CC = bwconncomp(vol);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        filtered_vol = false(size(vol));
        filtered_vol(CC.PixelIdxList{idx}) = true;
    
        % Find pixels that are above the threshold
        k = find(filtered_vol);
        [a,b,c]=ind2sub(size(densityGridSamp),k);
        edgeSize = max(size(densityGridSamp));
    
        keepers = 0*a; % Values that count as peaks
        for kk = 1:length(a)
            x1=a(kk);
            x2=b(kk);
            x3=c(kk);
            % This is an adaptation of the algorithm in pkfnd, adapted to 3D data
            if x1>1 && x1<edgeSize && x2>1 && x2<edgeSize && x3>1 && x3<edgeSize % avoid the image edges        
                if densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3) && ... % all 26 neighbors
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3) && ...    
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1-1,x2+1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2-1,x3+1) && ...
                   densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3-1) && densityGridSamp(x1,x2,x3) >= densityGridSamp(x1+1,x2+1,x3+1)
        
                        keepers(kk) = 1; % If the pixel is the brightest of all its neighbors, keep it
        
                end
        
            end
        end
        
        pkxSamp = [a(keepers>0) b(keepers>0) c(keepers>0)]; % Peaks in the mesh grid coordinate system
        
        %%% Calculate Fine Peaks Based on Point Cloud 
    
        sz = scale; %window for fitting neighbor peaks
    
        pkCSamp = NaN*ones(size(pkxSamp,1),4);
        for kk = 1:sum(keepers)
    
            r = sqrt(sum((dataSamp(:,1)-pkxSamp(kk,2)+meshPos).^2 + (dataSamp(:,2)-pkxSamp(kk,1)+meshPos).^2 + (dataSamp(:,3)-pkxSamp(kk,3)+meshPos).^2,2));
            neigh = (r<sz);
            if ~sum(neigh)
                error('No neighbors found around peak')
            end
    
            tmp = densitySamp(neigh);
            norm = sum(tmp);
            xavg = sum(dataSamp(neigh,1).*tmp)./norm; % Weighted average based on local density
            yavg = sum(dataSamp(neigh,2).*tmp)./norm;
            zavg = sum(dataSamp(neigh,3).*tmp)./norm;
    
            %  Radius of gyration (= sum(m_i * r_i^2)/sum(m_i)), m = mass, r = radius
            r = sqrt(sum((dataSamp(neigh,1)-xavg).^2 + (dataSamp(neigh,2)-yavg).^2 + (dataSamp(neigh,3)-zavg).^2,2));
            rg = sum(tmp.*r.^2)/norm;
    
            pkCSamp(kk,:) = [xavg yavg zavg rg]; % This is in nm relative in the superparticle coordinate system
    
        end
%         disp(pkCSamp)

%         figure(3)
%         set(gcf,'Position',[500 275 560*2 420*2])    
%         scatter3(xSamp,ySamp,zSamp, dSamp*5*10^6, dSamp, '.')
%         hold on
%         plot3(pkCSamp(:,1),pkCSamp(:,2),pkCSamp(:,3),'pr','MarkerFaceColor','r','MarkerSize',10);
%         hold off
%         colormap(cool)
%         xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
%         set(gca,'DataAspectRatio',[1 1 1])
%         set(gca,'FontSize',16)
%         colorbar       
%         set(gcf,'Color','white')
%         view(2)
%         drawnow
% 
%         dPeakSamp = mvksdensity(dataSamp,pkCSamp(:,1:3),'Bandwidth',scale/2); %  This is the probability density at the peaks
%         figure(4)
%         scatter3(dataSamp(:,1),dataSamp(:,2),densitySamp(:),[],densitySamp(:))
%         set(gca,'DataAspectRatio',[1 1 10^-6])
%         hold on
%         for kk = 1:size(pkCSamp,1)
%             plot3(pkCSamp(kk,1),pkCSamp(kk,2),dPeakSamp(kk),'pr','MarkerFaceColor','r','MarkerSize',10);
%         end
%         hold off
%         drawnow

        pkCLoo{ii} = pkCSamp;

        dataLoo{ii} = dataSamp;
        densityLoo{ii} = densitySamp;
        bladesLoo{ii} = bladesSamp;
        xLoo{ii} = xSamp;
        yLoo{ii} = ySamp;
        zLoo{ii} = zSamp;
        dLoo{ii} = dSamp;
        dPeakLoo{ii} = dPeakSamp;
        
%         toc

    end
    toc

    %% Measure Uncertainties
    loo = cell2mat(pkCLoo);
    looID = kmeans(looTest,size(pkC,1));
    looClus = NaN*ones(max(looID),3);
    for cc = 1:max(looID(:))
        rowsNow = looID == c;
        looClus(cc,:) = mean(loo(rowsNow,1:3));
    end
    allPeaks = cell2mat(pkCBoot);
    allH = cell2mat(dPeakBoot);
    idClus = dbscan(allPeaks(:,1:3),1,250); % epsilon = 1 means jumping the pixel-locking gaps, minpoints ~ proportion of instances, so 100 = ~10%

    figure(10)
    set(gcf,'Position',[400 400 560*2 420*2])
    scatter3(allPeaks(:,1),allPeaks(:,2),allPeaks(:,3),6,allH,'filled')
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
    h = colorbar;
    set(get(h,'Label'),'String','Peak Height')    
    saveas(gcf,[directory saveTag 'BootstrapAllPeaks.png'],'png')

    figure(11)    
    set(gcf,'Position',[800 300 560*2 420*2])
    scatter3(allPeaks(:,1),allPeaks(:,2),allPeaks(:,3),6,idClus,'filled')
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
    colormap([0.75 0.75 0.75 ; 0.75 0.75 0.75 ; lines(max(idClus(:)))])
    caxis([-1 max(idClus(:))])
    h = colorbar;
    set(get(h,'Label'),'String','DBSCAN cluster')    
    saveas(gcf,[directory saveTag 'BootstrapDBSCAN.png'],'png')
 %%   
    cClus = NaN*ones(max(idClus),3);
    alphaAll = NaN*ones(2,4,max(idClus));
    CIall = alphaAll;
    ahat = NaN*ones(max(idClus),4);
    z0 = ahat;
    clus2Orig = NaN*ones(max(idClus),1);
    clus2Loo = clus2Orig;    
    for cc = 1:max(idClus)
        rowsNow = idClus == cc;
        cClus(cc,:) = mean(allPeaks(rowsNow,1:3));

        %%% connect this cluster to the appropriate original peak and loo cluster
        idOrig = knnsearch(pkC(:,1:3),cClus(cc,:));
        idLoo = knnsearch(looClus,cClus(cc,:));
        clus2Orig(cc) = idOrig;
        clus2Loo(cc) = idLoo;

        sliceTest = (looID == idLoo);
        sliceTest = looTest(sliceTest,:);
        Ui = (size(sliceTest,1)-1)*(mean(sliceTest)-sliceTest);
        ahat(cc,:) = sum(Ui.^3)./(sum(Ui.^2).^3/2)/6;
    
        z0(cc,:) = norminv(sum(allPeaks(rowsNow,:)<pkC(idOrig,:))/sum(rowsNow));
        c = 0.025;
    
        alpha(1,:) = normcdf( z0(cc,:) + (z0(cc,:) + norminv(c)) / (1 - ahat(cc,:).*(z0(cc,:) + norminv(c))) );
        alpha(2,:) = normcdf( z0(cc,:) + (z0(cc,:) + norminv(1-c)) / (1 - ahat(cc,:).*(z0(cc,:) + norminv(1-c))) );
        alphaAll(:,:,cc) = alpha;
    
        CI = NaN*alpha;
        for dd = 1:size(alpha,2)
            CI(:,dd) = prctile(allPeaks(rowsNow,dd),alpha(:,dd)*100);
        end
        CIall(:,:,cc) = CI;
%         disp(cc)
%         disp(CI)

        figure(cc)
        set(gcf,'Position',[200 200 560*2.5 420])
        subplot(1,3,1)
        histogram(allPeaks(rowsNow,1));
        hold on
        plot(pkC(idOrig,1)*[1 1],[0 300],'Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(1,1)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(2,1)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(mean(allPeaks(rowsNow,1))*[1 1],[0 300],'LineWidth',2)
        hold off        
        title('x')

        subplot(1,3,2)
        histogram(allPeaks(rowsNow,2));
        hold on
        plot(pkC(idOrig,2)*[1 1],[0 300],'Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(1,2)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(2,2)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(mean(allPeaks(rowsNow,2))*[1 1],[0 300],'LineWidth',2)
        hold off
        title('y')

        subplot(1,3,3)
        histogram(allPeaks(rowsNow,3));
        hold on
        plot(pkC(idOrig,3)*[1 1],[0 300],'Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(1,3)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(CI(2,3)*[1 1],[0 300],'--','Color',[0.8500    0.3250    0.0980],'LineWidth',2)
        plot(mean(allPeaks(rowsNow,3))*[1 1],[0 300],'LineWidth',2)
        hold off
        title('z')

        saveas(gcf,[directory saveTag 'Bootstrap_cluster' num2str(cc) '.png'],'png')

    end
%%
    figure(12)
    set(gcf,'Position',[1200 200 560*2 420*2])
    plot3(pkC(:,1),pkC(:,2),pkC(:,3),'ok','MarkerFaceColor','k','MarkerSize',6)
    hold on
    for cc = 1:size(CIall,3)
        plot3(CIall(:,1,cc),pkC(clus2Orig(cc),2)*[1 1],pkC(clus2Orig(cc),3)*[1 1],'k','LineWidth',2)
        plot3(pkC(clus2Orig(cc),1)*[1 1],CIall(:,2,cc),pkC(clus2Orig(cc),3)*[1 1],'k','LineWidth',2)
        plot3(pkC(clus2Orig(cc),1)*[1 1],pkC(clus2Orig(cc),2)*[1 1],CIall(:,3,cc),'k','LineWidth',2)
    end    
    hold off
    set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)

%     % Add distances with propagated errors
%     function [d,ud] = distUncert(x1,x2,u1,u2)
%         % x1, x2 = coordinates
%         % u1, u2 = upper confidence interval
%     end
%     saveas(gcf,[directory saveTag 'BootstrapAllPeaks.png'],'png')

%% Bootstrap distance measurements
% What we really want is the distance, so that should be the derived
% quantity, not the peak locations...

% allPeaks = cell2mat(pkCBoot);
% allH = cell2mat(dPeakBoot);
% idClus = dbscan(allPeaks,1,100); % epsilon = 1 means jumping the pixel-locking gaps, minpoints ~ proportion of instances, so 100 = ~10%

eachPeak = [NaN*allPeaks NaN*allPeaks(:,1:2)];
counter = 1;

for cc = 1:length(pkCBoot)
    slice = pkCBoot{cc};
    eachPeak(counter:counter+size(slice,1)-1,:) = [slice cc*ones(size(slice,1),1) dPeakBoot{cc}];
    counter = counter+size(slice,1);
end

idE = dbscan(eachPeak(:,1:3),1,250);

figure(20)
set(gcf,'Position',[400 400 560*2 420*2])
scatter3(eachPeak(:,1),eachPeak(:,2),eachPeak(:,3),6,idE,'filled')
set(gca,'DataAspectRatio',[1 1 1],'FontSize',16)
colormap([0.75 0.75 0.75 ; 0.75 0.75 0.75 ; lines(max(idClus(:)))])
caxis([-1 max(idClus(:))])
h = colorbar;
set(get(h,'Label'),'String','DBSCAN cluster')   

distAll = NaN*ones(6,6,length(pkCBoot));
for cc = 1:length(pkCBoot)

%     fprintf([num2str(cc) ','])
%     if ~mod(cc,25)
%         fprintf('\n')
%     end

    slice = eachPeak(:,5) == cc;
    idSlice = idE(slice);
    slice = eachPeak(slice,:);

    % Remove the noise cluster from consideration
    noiseClus = (idSlice == -1);
    slice = slice(~noiseClus,:);
    idSlice = idSlice(~noiseClus);

%     if length(unique(idSlice)) ~= length(idSlice)
%         disp(cc)
%         error('Extra Peaks')
%     end

    toPdist = NaN*ones(6,3);
    for dd = 1:6
        sliceNow = idSlice == dd;
        if sum(sliceNow) > 1
            sliceNow = slice(sliceNow,1:3);
            warning(['Extra Peaks detected for instance ' num2str(cc)])
            toPdist(dd,:) = mean(sliceNow); % Is this the right behavior for a duplicate?
        elseif sum(sliceNow)
            toPdist(dd,:) = slice(sliceNow,1:3);
        end
    end

    distAll(:,:,cc) = squareform(pdist(toPdist));

end

figure(21)
for cc = 1:6
    for dd = 1:6
        subplot(6,6,6*(dd-1)+cc)
        histogram(distAll(cc,dd,:))
    end
end


    %% Save the data    
    save([directory saveTag '_PeakFindingAll.mat']) % This is likely to be huge....

end


%% Clean up
% close all
disp('Batch Complete')