function [r0p] = initial_momentum_guess()
G=7; %gamma factor
gammaBeta = sqrt(G^2-1);
dbs = 1.1;
Ip = [1.0, -2.0, 1.0];

rx_goal = pi/2;
ry_goal = 1.1;
fun=@(quads)rotation_cost_function(quads, rx_goal, ry_goal, gammaBeta, dbs);
options = optimoptions ('fsolve', 'Display', 'none','Algorithm', 'levenberg-marquardt');
options.MaxIter = 50000 ;
options.MaxFunEvals = 5000000 ;
options.FunctionTolerance = 1e-50;
options.OptimalityTolerance = 1e-50;
options.StepTolerance = 1e-25;
options = optimoptions(@lsqnonlin,options);
lb = [-6.0, -6.0, -6.0];
ub = [6.0, 6.0, 6.0];
[Is,~]=lsqnonlin(fun,Ip,lb,ub,options);
Is = round(Is,2);
[~, ~, ~, ~, rotationsx, rotationsy] = SER_decomposition(Is, gammaBeta, dbs);        
format = '\nguess currents %.4f_%.4f_%.4f: rotX = %.4f, rotY = %.4f. Goal: rx = %.4f, ry = %.4f';
fprintf(format, Is(1), Is(2), Is(3), rotationsx, rotationsy, rx_goal, ry_goal)

%initial beam parameters
qe=-1.6e-19;
beamdiv=1e-5; %beam divergence approximately
Spotsize=10e-4; %for gaussian dist
G=7;
dG=1e-6; %energy spread rms control parameter.
zsig=3e-12; %bunch length in seconds
Qbeam=2e6*qe*0;
nps=50000; %number of particles

s_beam = [1e-6 0 0 0; 0 2e-8 0 0; 0 0 3e-6 0; 0 0 0 0.5e-8];
centroid_beam = [0 0 0 0];
path_to_gpt='C:\Users\vgwr6\OneDrive\Desktop\UCLA\General_Particle_Tracer\bin\';
path_to_in_out_file=' D:\UCLA\PEGASUS\tomography\Generate_dataset\' ; %Must adapt this to your pc

%% Create beam file
    x1 = normrnd(0, 3e-4, [nps/2,1]);
    Px1 = normrnd(-2e-4, 1e-4, [nps/2,1]);
    gauss1 = [Px1, x1];
    x2 = normrnd(0, 3e-4, [nps/2,1]);
    Px2 = normrnd(2e-4, 1e-4, [nps/2,1]);
    gauss2 = [Px2, x2];
    gauss = vertcat(gauss1,gauss2);

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
    
     %% Output to txt to gdf
    fileID = fopen('input_distribution.txt','w');
    fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s\n','x','y','z','Bx','By','Bz','nmacro');
    for i = 1:nps
    fprintf(fileID,'%16.12e %16.12e %16.12e %16.12e %16.12e %16.12e %16.12e\n', x(i), y(i), z(i), Bx(i), By(i), Bz(i), nmacro(i));
    end
    fclose(fileID);
    system(strcat(path_to_gpt,'asci2gdf -o start.gdf input_distribution.txt'));
    
    %% Looping over quad settings        
    clear data A B 
    I4=Is(1);
    I5=Is(2);
    I6=Is(3);
    [I4, I5, I6]
    [ R , ~, ~ ,~] = transport_matrix_ff_focusing_green_var(-Is,gammaBeta,dbs);
    %% Run tomography_test.in using 'system()' function.
    system(strcat(path_to_gpt,'gpt -o results.gdf ',path_to_in_out_file,'tomography_test_fixed.in I4=',num2str(I4),' I5=',num2str(I5),' I6=',num2str(I6),' GPTLICENSE=1476385047'),'-echo');
    %% Create images
    data=loadgdf('results.gdf') ; 
    phsp1 = [data{1}.data.x data{1}.data.Bx data{1}.data.y data{1}.data.By];
    filename = strcat('start.bmp');
    imsize = [800,800];
    pixelcal = 25e-6;
    psf = 80e-6;
    A = create_image_from_partcoord(phsp1, imsize, pixelcal, psf, 1, filename);

    phsp2 = [data{2}.data.x data{2}.data.Bx data{2}.data.y data{2}.data.By];
    filename = strcat('D:\UCLA\PEGASUS\tomography\Generate_dataset\momentum_guess_analysis\target_',num2str(I4),'_',num2str(I5),'_',num2str(I6),'.bmp');
    B = create_image_from_partcoord(phsp2, imsize, pixelcal, psf, 1, filename);
%%
pixel = 25e-6;
npart_seed= 6e6;
centerx=size(A,1)/2;
centery=size(A,2)/2;
r0 = size(A,1)*pixel/2;
r2 = rand(npart_seed,1)*r0*r0;
phi = rand(npart_seed,1)*2*pi;
xpos=sqrt(r2).*cos(phi);
ypos=sqrt(r2).*sin(phi);
xposint =centerx+floor(xpos./pixel)+1;
yposint =centery+floor(ypos./pixel)+1; 

prob = rand(npart_seed,1);
idx = sub2ind(size(B),yposint,xposint);
xpos(B(idx)<prob)=[];
ypos(B(idx)<prob)=[];
xpos = downsample(xpos, 5);
ypos = downsample(ypos, 5);

rms_x0p = std(xpos) / R(1,2);
r0p = sqrt(2) * rms_x0p
end
