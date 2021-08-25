%%%WARNING, RUN THIS AFTER YOU HAVE RUN Tomgraphy_testfun%%%%
screen=1;
% quads=[1.1,-3.2,2.15];%quad currents
quads=[-I4,-I5,-I6];
G=7; %gamma factor
gammaBeta=sqrt(G^2-1); %gamma times beta
L=1.375; %distance from 3rd quad center to final screen.
position=0.055; %initial screen position
dbs = 1.1;


%%%triplet_focusing_example is the matrix propagator%%%%
% [ R , R_all, z_all ,pos] = triplet_focusing_example(quads,gammaBeta,L,position);

[ R , R_all, z_all ,pos] = transport_matrix_ff_focusing_green_var(quads,gammaBeta,dbs);

%%% need to pull data from GPT simulation to construct the initial sigma
%%% matrix.
x=data{screen}.data.x;
xp=data{screen}.data.Bx;
y=data{screen}.data.y;
yp=data{screen}.data.By;
phsp=zeros(length(x),4);
phsp(:,1)=x;
phsp(:,2)=xp;
phsp(:,3)=y;
phsp(:,4)=yp;
sigma0=cov(phsp);


%Initialize arrays to append C and S matrix elements at each z
Cpointsx=zeros(1,length(R_all));
Spointsx=zeros(1,length(R_all));
Cpointsy=zeros(1,length(R_all));
Spointsy=zeros(1,length(R_all));
%Initialize arrays to append Spotsizes at each z
stdx=zeros(1,length(R_all));
stdxp=zeros(1,length(R_all));
stdy=zeros(1,length(R_all));
stdyp=zeros(1,length(R_all));


for j=1:length(Spointsx)
Cpointsx(j)=R_all(1,1,j);
Spointsx(j)=R_all(1,2,j);
Cpointsy(j)=R_all(3,3,j);
Spointsy(j)=R_all(3,4,j);
S_aux = R_all(:,:,j)*sigma0*R_all(:,:,j)';
stdx(j) = sqrt(S_aux(1,1));
stdxp(j) = sqrt(S_aux(1,2));
stdy(j) = sqrt(S_aux(3,3));
stdyp(j) = sqrt(S_aux(3,4));
end


%% plot the results
figure
plot(z_all,Cpointsx,'LineWidth',2)
hold on
plot(z_all,Cpointsy,'LineWidth',2)
legend('C_x(z)','C_y(z)','Location','best')
xlabel('z(m)')
set(gca,'FontSize',15)


figure
plot(z_all,Spointsx,'LineWidth',2)
hold on
plot(z_all,Spointsy,'LineWidth',2)
legend('S_x(z)','S_y(z)','Location','best')
xlabel('z(m)')
set(gca,'FontSize',15)

%% compare GPT to the matrix code
figure
plot(z_all,stdx,'LineWidth',2)
hold on
plot(z_all,stdy,'LineWidth',2)
% hold on
% plot(stats{1}.data.position(2:end),stats{1}.data.stdx(2:end),'--','LineWidth',2)
% hold on
% plot(stats{1}.data.position(2:end),stats{1}.data.stdy(2:end),'--','LineWidth',2)
legend('\sigma_x(z)','\sigma_y(z)','\sigma_x(z) (GPT)','\sigma_y(z) (GPT)','Location','best')
xlabel('z(m)')
set(gca,'FontSize',15)

%%
fin = R * [x'; xp'; y'; yp'];
finalx = fin(1,:);
finalxp = fin(2,:);
finaly = fin(3,:);
finalyp = fin(4,:);
%%
figure(6)
ax1=subplot(2,3,1)
scatter(finalx*1e+6,finalxp*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('\beta_x (mrad)')
set(gca,'FontSize',15)

ax2=subplot(2,3,2)
scatter(finaly*1e+6,finalyp*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('y (um)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)

ax3=subplot(2,3,3)
scatter(finaly*1e+6,finalxp*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('y (um)')
ylabel('\beta_x (mrad)')
set(gca,'FontSize',15)

ax5=subplot(2,3,4)
scatter(finalx*1e+6,finalyp*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)

ax5=subplot(2,3,5)
scatter(finalxp*1e+3,finalyp*1e+3,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('\beta_x (mrad)')
ylabel('\beta_y (mrad)')
set(gca,'FontSize',15)
axis equal

ax6=subplot(2,3,6)
scatter(finalx*1e+6,finaly*1e+6,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('y (um)')
set(gca,'FontSize',15)
axis equal

figure(7)
scatter(finalx*1e+6,finaly*1e+6,0.5)
ax = gca;
ax.YDir = 'normal'
xlabel('x (um)')
ylabel('y (um)')
set(gca,'FontSize',15)
axis equal

%%
fileloc='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\'
resfactor=1;
saveImage=0;
pixcal=13.614e-6; %13.6143-6 old calibration? not really sure where it came from.
if saveImage==1
     yag_Image = zeros(floor(800/resfactor),floor(800/resfactor));
     for l=1:length(finalx)
         xpos=floor(finalx(l)/(resfactor*pixcal)+400/resfactor);
         ypos=floor(-finaly(l)/(resfactor*pixcal)+400/resfactor);
         if xpos>0 && xpos<=800/resfactor && ypos>0 && ypos<=800/resfactor
             yag_Image(ypos,xpos)=yag_Image(ypos,xpos)+1;
         end 
     end
     %renormalize and smooth image
     maxIm=max(max(yag_Image));
     yag_Image=yag_Image.*floor(400/maxIm);
     radius=1; %pixels
     filter=fspecial('disk',radius);
     yag_Image=conv2(yag_Image,filter);
     yag_Image = imresize(yag_Image,[800 800]);
     yag_Image=uint8(yag_Image);
     imageName = strcat('matrix_image_1.1_-3.2_2.15','.bmp');
     imwrite(yag_Image,strcat(fileloc,imageName));
end
final_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\matrix_image_1.1_-3.2_2.15.bmp');
% sum(final_image,'all')
% final_image = (100000 / sum(final_image,'all')) * final_image;
figure
imagesc(yag_Image)
colorbar
