clear all
close all
%%
%initial beam parameters
qe=-1.6e-19;
beamdiv=1e-5; %beam divergence approximately
Spotsize=10e-4; %for gaussian dist
G=7;
dG=1e-6; %energy spread rms control parameter.
zsig=3e-12; %bunch length in seconds
Qbeam=2e6*qe*0;
nps=100000; %number of particles

s_beam = [1e-6 0 0 0; 0 2e-8 0 0; 0 0 3e-6 0; 0 0 0 0.5e-8];
centroid_beam = [0 0 0 0];
%path_to_gpt='C:\Users\vgwr6\OneDrive\Desktop\UCLA\General_Particle_Tracer\bin\';
path_to_in_out_file=pwd; %Must adapt this to your pc

%% 
x1 = normrnd(0, 3e-4, [nps/2,1]);
Px1 = normrnd(-4e-4, 1e-4, [nps/2,1]);
gauss1 = [Px1, x1];
x2 = normrnd(0, 3e-4, [nps/2,1]);
Px2 = normrnd(4e-4, 1e-4, [nps/2,1]);
gauss2 = [Px2, x2];
gauss = vertcat(gauss1,gauss2);
% gauss = gauss/norm(gauss);
%% Create beam file
    nmacro  = Qbeam/qe/nps * ones(nps,1);
    TPS     = mvnrnd(centroid_beam,s_beam,nps);
    x       = TPS(:,1);
%     Bx      = TPS(:,2);
    Bx      = gauss(:,1);
    y       = TPS(:,3);
%     By      = TPS(:,4);
    By      = gauss(:,2);
    z       = normrnd(3.1,zsig, [nps 1]) - 3*zsig;
    G       = normrnd(G , dG , [nps 1]);
    B       = sqrt(1-1./G.^2);
    Bz      = sqrt(B.^2-Bx.^2-By.^2);
    
    Bx = Bx - x/5;
    By = By + y/3;
    

    %% Output to txt to gdf
    fileID = fopen('input_distribution.txt','w');
    fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s\n','x','y','z','Bx','By','Bz','nmacro');
    for i = 1:nps
    fprintf(fileID,'%16.12e %16.12e %16.12e %16.12e %16.12e %16.12e %16.12e\n', x(i), y(i), z(i), Bx(i), By(i), Bz(i), nmacro(i));
    end
    fclose(fileID);
    system(strcat('asci2gdf -o start.gdf input_distribution.txt'));
    
%% Quad settings (currents) 

imsize = [300,300];
pixelcal = 40e-6;
psf = 40e-6;
current_data = importdata(strcat(pwd,'\..\Experiment\current_settings.csv'));
current_fin = current_data.data;
quadlist(:,1) = current_fin(:,1);
quadlist(:,2) = current_fin(:,2);
quadlist(:,3) = current_fin(:,3);
%% Looping over quad settings        
for iq = 1:size(quadlist)
    clear data A B 
    I4=quadlist(iq,1);
    I5=quadlist(iq,2);
    I6=quadlist(iq,3);
    fprintf(1,'settings %i out of %i : quadvalues %3.3f %3.3f %3.3f \n',iq, size(quadlist,1), I4, I5 ,I6);
    %% Run tomography_test.in using 'system()' function.
    system(strcat('gpt -o results.gdf tomography_test_struct.in I4=',num2str(I4),' I5=',num2str(I5),' I6=',num2str(I6),' GPTLICENSE=1476385047'),'-echo');
    %% Create images
    data=loadgdf('results.gdf') ; 
    if (iq == 1)     
        phsp1 = [data{1}.data.x data{1}.data.Bx data{1}.data.y data{1}.data.By];
        figure()
        scatter(data{1}.data.Bx, data{1}.data.By, 5, 'MarkerFaceAlpha', 0.05);
        title('Initial momentum distribution');
        filename = strcat('start_dg.bmp');
        A = create_image_from_partcoord(phsp1, imsize, pixelcal, psf, 1, filename);
    end
        phsp2 = [data{2}.data.x data{2}.data.Bx data{2}.data.y data{2}.data.By];
        filename = strcat(pwd,'\dg_gpt\target_',num2str(I4),'_',num2str(I5),'_',num2str(I6),'.bmp');
        B = create_image_from_partcoord(phsp2, imsize, pixelcal, psf, 1, filename);
end

