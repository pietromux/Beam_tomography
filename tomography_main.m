clear all
close all

%initial beam parameters
qe=-1.6e-19;
beamdiv=1e-5; %beam divergence approximately
Spotsize=10e-4; %for gaussian dist
dG=1e-6; %energy spread rms control parameter.
ztime=10e-12; %bunch length in seconds
G=7;
Qbeam=2e+6*qe*0;
numstd=0.25;
npar=50000; %numbber of particles

%Quad settings (currents) arbitrary for now%
I3=1.1;
I4=-3.2;
I5=2.15;
I6=5*0;

%what are the phase advances for arbitrary settings?

%% %%%%%%%%%%%%%%%%%%%%%%%%%Run GPT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
infile = 'Tomography_test.in';
gdffile = 'Tomography_test_solution.gdf';
path_to_gpt='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\General_Particle_Tracer\bin\'; %Must adapt this to your pc
path_to_in_out_file=' C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\' ; %Must adapt this to your pc
license=' GPTlicense=1476385047';

callGPT = strcat(path_to_gpt,'gpt -v -o ',path_to_in_out_file,gdffile,path_to_in_out_file,infile,...
' Qtot=',num2str(Qbeam),' npar=',num2str(npar),' Spotsize=',num2str(Spotsize),' beamdiv=',num2str(beamdiv),' G=',num2str(G),' dG=',num2str(dG),...
' ztime=',' I3=',num2str(I3),' I4=',num2str(I4),' I5=',num2str(I5),' I6=',num2str(I6),num2str(ztime),license)

system(callGPT,'-echo') ;
data=load_gdf(strcat(path_to_in_out_file,gdffile)) ; 
callGPT2=strcat(path_to_gpt,'gdfa -o',path_to_in_out_file,'result.gdf ',path_to_in_out_file,gdffile,' position ',' stdBy ',' stdBx ',' stdBz',' stdx ',' stdy',' CSgammax',' CSgammay',' CSalphax',' CSalphay', ' CSbetax',' CSbetay',' nemixrms',' nemiyrms',' nemirrms',' nemizrms',' avgz',' stdz',' avgG',' stdG');
system(callGPT2,'-echo') ;
stats=load_gdf(strcat(path_to_in_out_file,'result.gdf')); 

%% Initial Yag Screen
u=zeros(2,length(data(1).d.x));
x2ps=data(1).d.x;
y2ps=data(1).d.y;

fileloc='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\'
resfactor=1;
saveImage=1;
pixcal=13.614e-6; % 13.614e-6 old calibration? not really sure where it came from.

if saveImage==1
     yag_Image_init = zeros(floor(800/resfactor),floor(800/resfactor));
     for l=1:length(x2ps)
         xpos=floor(x2ps(l)/(resfactor*pixcal)+400/resfactor);
         ypos=floor(-y2ps(l)/(resfactor*pixcal)+400/resfactor);
         if xpos>0 && xpos<=800/resfactor && ypos>0 && ypos<=800/resfactor
             yag_Image_init(ypos,xpos)=yag_Image_init(ypos,xpos)+1;
         end
     end
     
     %renormalize and smooth image
     maxIm=max(max(yag_Image_init));
     yag_Image_init=yag_Image_init.*floor(400/maxIm);
     radius=1; %pixels
     filter=fspecial('disk',radius);
     yag_Image_init = conv2(yag_Image_init,filter);
     yag_Image_init = imresize(yag_Image_init,[800 800]);
     yag_Image_init = uint8(yag_Image_init);
     imageName = strcat('GPT_image_init2_1.1_-3.2_2.15','.bmp');
     imwrite(yag_Image_init,strcat(fileloc,imageName));
end



init_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\GPT_image_init_1.1_-3.2_2.15.bmp');
numpart = sum(init_gpt_image, 'all');
figure()
imagesc(init_gpt_image)
colorbar




%% Final Yag Screen
u=zeros(2,length(data(end).d.x));
x2ps=data(end).d.x;
y2ps=data(end).d.y;

fileloc='C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\'
resfactor=1;
saveImage=1;
pixcal=13.614e-6; %old calibration? not really sure where it came from.
if saveImage==1
     yag_Image = zeros(floor(800/resfactor),floor(800/resfactor));
     for l=1:length(x2ps)
         xpos=floor(x2ps(l)/(resfactor*pixcal)+400/resfactor);
         ypos=floor(-y2ps(l)/(resfactor*pixcal)+400/resfactor);
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
     imageName = strcat('GPT_image_init2_1.1_-3.2_2.15','.bmp');
     imwrite(yag_Image,strcat(fileloc,imageName));
end

final_gpt_image=imread('C:\Users\vgwr6\Desktop\UCLA\soph\Musumeci_Lab\GPT\Test_Images\GPT_image_init2_1.1_-3.2_2.15.bmp');
final_gpt_image = 1/sum(final_gpt_image, 'all')*double(final_gpt_image);
figure()
imagesc(final_gpt_image)
colorbar




%% first iteration
[init4DCoordOld, fin4DCoord, finIm, finImDisp, diffMat, posCoord, negCoord, sigmaRecon, stReconMeasXY] = first_propagation(init_gpt_image, final_gpt_image, R, numpart, 3e-20, pixcal, resfactor);
sum(diffMat, 'all')
%% algorithm loop
for i = 1 : 5
    init4DCoord = backward_propagation(init4DCoordOld, fin4DCoord, diffMat, pixcal, resfactor);
    [init4DCoordOld, fin4DCoord, finIm, finImDisp, diffMat] = forward_propagation(init4DCoord, final_gpt_image, R, pixcal, resfactor);
    figure()
    imagesc(finImDisp)
    colorbar
    figure()
    imagesc(diffMat)
    colorbar
    ConvHist = []
    ConvHist(end+1) = sum(diffMat, 'all')
    i
end