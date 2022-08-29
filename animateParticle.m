function animateParticle(data,scale,figName,directory,axLim,cLim)

dAng = 5;
tPause = 10;

density = mvksdensity(data,data,'Bandwidth',scale/2); %  This is the probability density

figure(1)
set(gcf,'Position',[500 275 560*2 420*2])    
scatter3(data(:,1),data(:,2), data(:,3), density*5*10^6, density, '.')
colormap(cool)
xlabel('x (nm)'),ylabel('y (nm)'),zlabel('z (nm)')
set(gca,'DataAspectRatio',[1 1 1])
set(gca,'FontSize',16)
colorbar       
set(gcf,'Color','white')

if exist('axLim','var') && ~isempty(axLim)
    axis(axLim)
    addTag = 'manualAxis';
else
    axMin = floor(min(data(:)));
    axMax = ceil(max(data(:)));
    axis(repmat([axMin axMax],[1 3]))
    addTag = '';
end
if exist('cLim','var') && ~isempty(cLim)
    caxis(cLim)
end

saveas(gcf,[directory filesep figName '_' addTag 'coolColormapPDF.png'],'png')  

%%%% Rotate and record
v = VideoWriter([directory figName '_pdf' addTag '.mp4'],'MPEG-4');
v.FrameRate = 15;
open(v)

% f = getframe(gcf);
% writeVideo(v,f)

axis vis3d
for ii = 1:360/dAng % rotate in theta (horizontal)
    camorbit(dAng,0)
    f = getframe(gcf);
    writeVideo(v,f);
end
for ii = 1:tPause % pause before switching axes
    writeVideo(v,f);
end
for ii = 1:360/dAng % rotate i phi (veritical)
    camorbit(0,dAng)
    f = getframe(gcf);
    writeVideo(v,f);
end
for ii = 1:tPause % pause before ending
    writeVideo(v,f);
end

close (v)
