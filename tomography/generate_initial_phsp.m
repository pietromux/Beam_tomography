% initial parameters
clear phsp
clear msp

%%
dG=1e-6; %energy spread rms control parameter.
ztime=10e-12; %bunch length in seconds


npart_seed= 6e6; % number of seed particles to sample initial distribution

if (sim)
    filepath = 'D:\UCLA\PEGASUS\tomography\';
    filename = strcat(filepath,'Generate_dataset\start.bmp');
    pixel = 40e-6;
    G=7;
else
    filepath = 'D:\UCLA\PEGASUS\8_19_21_crop\';
    filename = strcat(filepath, 'start.bmp');
    pixel = 10.59322e-6;
    G=17;
end
A = imread(filename);
figure(202)
imagesc(A)
colorbar
title('Initial spatial distribution')

centerx=size(A,2)/2;
centery=size(A,1)/2;
r0x = size(A,2)*pixel/2;
r0y = size(A,1)*pixel/2;

% create distribution function in the interval (0,1) from image profile
B = double(A);
peak=max(max(B));
B = B/peak;
B = B;
%% look into this part a bit more and understand
% Extract particle coordinates with initscreen spatial distribution
% r2 = rand(npart_seed,1)*r0*r0;
% phi = rand(npart_seed,1)*2*pi;
% xpos=sqrt(r2).*cos(phi);
% ypos=sqrt(r2).*sin(phi);
xpos = rand(npart_seed,1) * r0x * 2 - r0x;
ypos = rand(npart_seed,1) * r0y * 2 - r0y;
xposint =floor(centerx+xpos./pixel)+1;
yposint =floor(centery-ypos./pixel)+1; 

prob = rand(npart_seed,1);
idx = sub2ind(size(B),yposint,xposint);
xpos(B(idx)<prob)=[];
ypos(B(idx)<prob)=[];
dratio = round(length(xpos)/50500);
xpos = downsample(xpos, dratio);
ypos = downsample(ypos, dratio);
phsp(:,1)=xpos';
phsp(:,3)=ypos';

npart = size(xpos,1)
%% Set random angular distribution
if (strcmp(flagdist,'gaussian'))    
    phsp(:,2) = r0p*randn(npart,1);
    phsp(:,4) = r0p*randn(npart,1);
elseif (strcmp(flagdist,'good guess'))
    phsp(:,2) = 4e-4*randn(npart,1);
    phsp(:,4) = 2e-4*randn(npart,1);
else   
    r2p = rand(npart,1)*r0p*r0p;
    phip = rand(npart,1)*2*pi;
    phsp(:,2)=sqrt(r2p).*cos(phip);
    phsp(:,4)=sqrt(r2p).*sin(phip);
end

phsp(:,5) = randn(npart,1)*ztime;
energy = randn(npart,1)*dG+G;
phsp(:,6) = sqrt(energy.^2 - 1);

%% Visualize initial phase space
figure(300)
sgtitle('Initial Guess')
subplot(3,2,1)
dscatter(phsp(:,1),phsp(:,3))
xlim([1.1*min(min(phsp(:,1))) 1.1*max(max(phsp(:,1)))]);
ylim([1.1*min(min(phsp(:,3))) 1.1*max(max(phsp(:,3)))]);
xlabel('x[m]')
ylabel('y[m]')
colorbar
subplot(3,2,2)
dscatter(phsp(:,2),phsp(:,4))
xlim([1.1*min(min(phsp(:,2))) 1.1*max(max(phsp(:,2)))]);
ylim([1.1*min(min(phsp(:,4))) 1.1*max(max(phsp(:,4)))]);
xlabel('xp[rad]')
ylabel('yp[rad]')
colorbar
subplot(3,2,3)
dscatter(phsp(:,1),phsp(:,2))
xlim([1.1*min(min(phsp(:,1))) 1.1*max(max(phsp(:,1)))]);
ylim([1.1*min(min(phsp(:,2))) 1.1*max(max(phsp(:,2)))]);
xlabel('x[m]')
ylabel('xp[rad]')
colorbar
subplot(3,2,4)
dscatter(phsp(:,3),phsp(:,4))
xlim([1.1*min(min(phsp(:,3))) 1.1*max(max(phsp(:,3)))]);
ylim([1.1*min(min(phsp(:,4))) 1.1*max(max(phsp(:,4)))]);
xlabel('y[m]')
ylabel('yp[rad]')
colorbar
subplot(3,2,5)
dscatter(phsp(:,1),phsp(:,4))
xlim([1.1*min(min(phsp(:,1))) 1.1*max(max(phsp(:,1)))]);
ylim([1.1*min(min(phsp(:,4))) 1.1*max(max(phsp(:,4)))]);
xlabel('x[m]')
ylabel('yp[rad]')
colorbar
subplot(3,2,6)
dscatter(phsp(:,3),phsp(:,2))
xlim([1.1*min(min(phsp(:,3))) 1.1*max(max(phsp(:,3)))]);
ylim([1.1*min(min(phsp(:,2))) 1.1*max(max(phsp(:,2)))]);
xlabel('y[m]')
ylabel('xp[rad]')
colorbar
sigma_guess=cov(phsp(:,1:4))
