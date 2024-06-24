%% Prep and Running CBathy
clear all
close all

%% Directory of Rectified Frames
% fdir='/Users/dylananderson/Documents/projects/satelib/s108stabilize.avi'; %Your output from A_export_frames 
% fdir='/Users/dylananderson/Documents/projects/satelib/s108_20230511T190757Z/stabilizedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s115_20230419T191912Z/unstretchedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s111_20230605T191503Z/unstretchedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s3_20230625T154019Z/stretchedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s112_20230905T145032Z/stretchedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s112_20230802T145627Z/unstretchedFrames/'; 
%fdir='/Users/dylananderson/Documents/projects/satelib/s108_20230511T190757Z/unstretchedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s112_20230905T145032Z/unstretchedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s107_20230919T191836Z/unstretchedFrames/'; 
% fdir='/Users/dylananderson/Documents/projects/satelib/s112_20230802T145627Z/unstretchedFrames/';
fdir='/Users/jonathan/Desktop/frf_collects/s111_20230605T191503Z/frames/';

cd(fdir)
dc = dir('*.tiff'); % loads all the image infos

Xv = -700:1:1299;
Yv = -300:1:1699;
[X,Y]=meshgrid(Xv,Yv);

% tRes = 15;
% t = (1/2):(1/2):55;
% tRes = 10;
% t = (1/3):(1/3):60;
% tRes = 6;
% t = (1/5):(1/5):30;
% tRes = 5;
t = (1/6):(1/6):60;

% tRes = 3;
% t = 1/10:(1/10):30;

% tRes = 2;
% t = 1/15:(1/15):60;

tRes = 1;
t = 1/30:(1/30):30;

% Grid for s111_20230605T191503Z
res = 2;
sx = 700;
ex = 1800;
sy = 300;
ey = 1300;
[Xg,Yg] = meshgrid(Xv(sx:res:ex),Yv(sy:res:ey));

% % Grid for s111_20230419
% res = 2;
% sx = 700;
% ex = 1800;
% sy = 800;
% ey = 2000;
% [Xg,Yg] = meshgrid(Xv(sx:res:ex),Yv(sy:res:ey));

% %Display As Example
I=flipud(imread(fullfile(fdir,dc(1).name)));
%Ig=uint8(interp2(X,Y,double(I),Xg,Yg));
Ig = I(sy:res:ey,sx:res:ex);

% figure
% pcolor(X,Y,I(:, :, 1)) % needed to modify this to be only one dim, is I supposed to be 3D? - JCM
% shading flat
% xlabel('X [m]')
% ylabel('Y [m]')
% colormap(gray)
% hold on
% % rectangle('Position',[min(min(Xg)) min(min(Yg))  max(max(Xg))-min(min(Xg)) max(max(Yg))-min(min(Yg)) ],'edgecolor','g')
% close all
% 
TS = Ig;
% c = 1;
% for i = 1:tRes:1800%length(dc) % loop through all your images to resize
%     % load image
%     baseFileName = dc(i).name;
%     image = flipud(imread(baseFileName));
%     %resize image
%     TS(:,:,c) = adapthisteq(image(sy:res:ey,sx:res:ex));
%     pcolor(Xg,Yg,TS(:,:,c))
%     shading flat
%     xlabel('X [m]')
%     ylabel('Y [m]')
%     colormap(gray)
%     axis equal
%     title(num2str(c))
%     pause(1/900)
%     c = c + 1;
% end

%%
% Cut T to indices used above
tepoch_cut=t;%(t-datenum(1970,1,1)).*24*3600; %cbathy likes epoch time


%Display As Example
figure
pcolor(Xg,Yg,TS)
shading flat
xlabel('X [m]')
ylabel('Y [m]')
colormap(gray)
hold on
rectangle('Position',[min(min(Xg)) min(min(Yg))  max(max(Xg))-min(min(Xg)) max(max(Yg))-min(min(Yg)) ],'edgecolor','g')


%% Get Into CBathy Format
% This section of code prepares data for use with the CBathy bathymetry estimation tool. It assumes all grid points have valid depth measurements and reshapes the data to fit CBathy's input format.

% Format Size
[r, c] = size(TS(:, :, 1));  % Get number of rows (r) and columns (c) from the first time slice of TS

% Create Grid Indices
[xindgrid, yindgrid] = meshgrid(1:c, 1:r);  % Create 2D grids of x and y indices for all points in the grid

% Convert to Linear Indices
rowIND = yindgrid(:);   % Flatten row indices into a single column vector
colIND = xindgrid(:);   % Flatten column indices into a single column vector

% Reshape Data
counter = 0;
data = [];  % Initialize empty data matrix
for i = 1:length(rowIND(:))  % Loop through each grid point
    counter = counter + 1;
    % Extract time series for this point, reshape into column vector, and append to data
    data(:, counter) = reshape(TS(rowIND(i), colIND(i), :), length(tepoch_cut), 1); 
end

%% Prepare Data for CBathy

% XYZ Coordinates (Assuming 0 for Z to exclude tides)
xyz = [Xg(:) Yg(:) Xg(:) * 0];  % Create matrix of XYZ coordinates. Xg and Yg are assumed to be vectors of X and Y coordinates.

% Load CBathy Configuration (from 'argus02b.m' file with modifications)
argus02b;  
bathy.params = params;  % Apply configuration to bathymetry parameters

% Set Geographic Boundaries for the Region of Interest
bathy.params.xyMinMax = [0, 900, 0, 1000];  % Define min/max X and Y coordinates

% Camera Information (Required for CBathy, even with just one 'camera')
cam = ones(size(xyz, 1), 1);  % Create a dummy camera matrix (all values are 1)

% Optional: Filter Frequency Bands (If needed)
% bathy.params.fB = bathy.params.fB(4:10);  

%% Run Cbaty- From Latest Master Branch
%bathy = analyzeBathyCollect(xyz, tepoch_cut, double(data), cam, bathy)
[f, G, bathy] = prepBathyInputShort(xyz, tepoch_cut, double(data), bathy);

bathy = analyzeBathyCollectShort(xyz, tepoch_cut, double(data), cam, bathy);

close all

plotStacksAndPhaseMaps(xyz,t,data,f,G, params)


%% Plot Cbathy
f1=figure;
% 
% subplot(121)
% pcolor(Xg,Yg,TS(:,:,1))
% colormap("gray")
% % pcolor(bathy.xm,bathy.ym,bathy.timex)
% shading flat
% set(gca,'ydir','normal')
% hold on
% shading flat
% xlabel('X [m]')
% ylabel('Y [m]')
% 

a=gca; 
% pcolor(X(1,:),Y(:,1),imread(fullfile(fdir,dc(1).name)))
pcolor(Xg,Yg,TS(:,:,1))
% pcolor(bathy.xm,bathy.ym,bathy.timex)
colormap(a,"gray")
shading flat
set(gca,'ydir','normal')
hold on
a2=axes;
linkaxes([a a2],'xy')
a2.Color='none'
hold on
pcolor(bathy.xm,bathy.ym,-bathy.fCombined.h)
shading flat
h=colorbar('south')
caxis([-10 0])
h.Color='w'
h.Label.String='Depth [m]';
xlabel('X [m]')
ylabel('Y [m]')


figure(13)
% subplot(122)

i = 4;
ind = find(abs(f-params.fB(i)) == min(abs(f-params.fB(i))));
% subplot(nRows, nCols, i,'FontSize',7); hold on
h=scatter(xyz(:,1),xyz(:,2),3,angle(G(ind,:)),'filled');
xlabel('x (m)'); 
ylabel('y (m)'); 
axis equal;
caxis([-pi pi]);
axis ([ min(xyz(:,1)) max(xyz(:,1)) min(xyz(:,2)) max(xyz(:,2))]);
view(2); 
title(['f = ' num2str(params.fB(i),'%0.3g') ' Hz'],'FontWeight','normal','FontSize',9); 
grid on


save("s112_20230802_cBathy2Hz60s.mat","bathy","data","params","f","G","Xg","Yg","TS","xyz","cam","tepoch_cut",'-v7.3')


% 
% figure(14)
% subplot(211)
% plot(xFRF,elevation(:,13),'linewidth',2)
% hold on
% plot(bathy.xm(4:end),-bathy.fCombined.h(21,4:end),'linewidth',2)
% ylim([-9, 0])
% title('yFRF = 150')
% ylabel('NAVD88 (m)')
% % legend('Survey June 7th','cBathy 3.0')
% xlim([100, 1000])
% subplot(212)
% plot(xFRF,elevation(:,24),'linewidth',2)
% hold on
% plot(bathy.xm(3:end),-bathy.fCombined.h(46,3:end),'linewidth',2)
% ylim([-9, 0])
% legend('Survey June 7th','cBathy 3.0 June 5th')
% title('yFRF = 450')
% ylabel('NAVD88 (m)')
% xlabel('xFRF (m)')