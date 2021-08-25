% initial parameters
clear phsp
clear msp

%%
dG=2e-6; %energy spread rms control parameter.
ztime=10e-12; %bunch length in seconds
G=7.807;
pixel = 40e-6;
npart_seed= 1e6; % number of seed particles to sample initial distribution

flagdist = 'uniform';
filepath = pwd;
filename = strcat(filepath,'\Generate_dataset\start_dg.bmp');
A = imread(filename);
figure()
imagesc(A)
colorbar
title('Initial spatial distribution')

centerx=size(A,1)/2;
centery=size(A,2)/2;
r0 = size(A,1)*pixel/2;

% create distribution function in the interval (0,1) from image profile
B = double(A);
peak=max(max(B));
B = B/peak;
B = B;
%% Generate particles in random squadre and removes particles to match B image
% Extract particle coordinates with initscreen spatial distribution
% r2 = rand(npart_seed,1)*r0*r0;
% phi = rand(npart_seed,1)*2*pi;
% xpos=sqrt(r2).*cos(phi);
% ypos=sqrt(r2).*sin(phi);
xpos = rand(npart_seed,1) * r0 * 2 - r0;
ypos = rand(npart_seed,1) * r0 * 2 - r0;
xposint =centerx+floor(xpos./pixel)+1;
yposint =centery+floor(-ypos./pixel)+1; 

prob = rand(npart_seed,1);
idx = sub2ind(size(B),yposint,xposint);
xpos(B(idx)<prob)=[];
ypos(B(idx)<prob)=[];
%xpos = downsample(xpos, 5);
%ypos = downsample(ypos, 5);
phsp(:,1)=xpos';
phsp(:,3)=ypos';
npart = size(xpos,1)

%% Set random angular distribution
if (strcmp(flagdist,'gaussian'))    
    phsp(:,2) = r0p*randn(npart,1);
    phsp(:,4) = r0p*randn(npart,1);
else   
    r2p = rand(npart,1)*r0p*r0p;
    phip = rand(npart,1)*2*pi;
    phsp(:,2)=sqrt(r2p).*cos(phip);
    phsp(:,4)=sqrt(r2p).*sin(phip);
end

phsp(:,5) = randn(npart,1)*ztime;
energy = randn(npart,1)*dG+G;
phsp(:,6) = sqrt(energy.^2 - 1);

%% Visualize initial phase space guess
figure()
sgtitle('Initial Guess')
subplot(2,2,1)
dscatter(phsp(:,1),phsp(:,3))
xlim([1.1*min(min(phsp(:,1))) 1.1*max(max(phsp(:,1)))]);
ylim([1.1*min(min(phsp(:,3))) 1.1*max(max(phsp(:,3)))]);
title('x vs. y')
colorbar
subplot(2,2,2)
dscatter(phsp(:,2),phsp(:,4))
xlim([1.1*min(min(phsp(:,2))) 1.1*max(max(phsp(:,2)))]);
ylim([1.1*min(min(phsp(:,4))) 1.1*max(max(phsp(:,4)))]);
title('xp vs. yp')
colorbar
subplot(2,2,3)
dscatter(phsp(:,1),phsp(:,2))
xlim([1.1*min(min(phsp(:,1))) 1.1*max(max(phsp(:,1)))]);
ylim([1.1*min(min(phsp(:,2))) 1.1*max(max(phsp(:,2)))]);
title('x vs. xp')
colorbar
subplot(2,2,4)
dscatter(phsp(:,3),phsp(:,4))
xlim([1.1*min(min(phsp(:,3))) 1.1*max(max(phsp(:,3)))]);
ylim([1.1*min(min(phsp(:,4))) 1.1*max(max(phsp(:,4)))]);
title('y vs. yp')
colorbar

%% Visualize actual initial phase space

start = loadgdf(strcat(filepath,'\Generate_dataset\struct.gdf'));
phsp_start = [start{1}.data.x start{1}.data.Bx start{1}.data.y start{1}.data.By];
s_beam_target = cov(phsp_start);

figure()
sgtitle('Actual Initial Distribution')
subplot(2,2,1)
dscatter(phsp_start(:,1),phsp_start(:,3))
xlim([1.1*min(min(phsp(:,1))) 1.1*max(max(phsp(:,1)))]);
ylim([1.1*min(min(phsp(:,3))) 1.1*max(max(phsp(:,3)))]);
title('x vs. y')
colorbar
subplot(2,2,2)
dscatter(phsp_start(:,2),phsp_start(:,4))
xlim([1.1*min(min(phsp(:,2))) 1.1*max(max(phsp(:,2)))]);
ylim([1.1*min(min(phsp(:,4))) 1.1*max(max(phsp(:,4)))]);
title('xp vs. yp')
colorbar
subplot(2,2,3)
dscatter(phsp_start(:,1),phsp_start(:,2))
xlim([1.1*min(min(phsp(:,1))) 1.1*max(max(phsp(:,1)))]);
ylim([1.1*min(min(phsp(:,2))) 1.1*max(max(phsp(:,2)))]);
title('x vs. xp')
colorbar
subplot(2,2,4)
dscatter(phsp_start(:,3),phsp_start(:,4))
xlim([1.1*min(min(phsp(:,3))) 1.1*max(max(phsp(:,3)))]);
ylim([1.1*min(min(phsp(:,4))) 1.1*max(max(phsp(:,4)))]);
title('y vs. yp')
colorbar

