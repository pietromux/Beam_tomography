clear all
close all
%%
current_settings = readtable('D:\UCLA\PEGASUS\tomography\current_settings_G7_no_constraint.xls');
current_settings_cell_array = current_settings;
current_settings_cell_array(1,:) = [];
current_settings_matrix = table2array(current_settings_cell_array);
%%
current_settings_matrix(:,8) = current_settings_matrix(:,6) - current_settings_matrix(:,4);
current_settings_matrix(:,9) = current_settings_matrix(:,7) - current_settings_matrix(:,5);

diff_rx = current_settings_matrix(:,8);
diff_ry = current_settings_matrix(:,9);

I4_array = current_settings_matrix(:,1);
I5_array = current_settings_matrix(:,2);
I6_array = current_settings_matrix(:,3);

test_settings(:,1) = I4_array(diff_rx<1e-3 & diff_ry<1e-3 & abs(I4_array)<6);
test_settings(:,2) = I5_array(diff_rx<1e-3 & diff_ry<1e-3 & abs(I4_array)<6);
test_settings(:,3) = I6_array(diff_rx<1e-3 & diff_ry<1e-3 & abs(I4_array)<6);
test_settings = round(test_settings, 2);
test_settings = unique(test_settings, 'rows');
%%
qe=-1.6e-19;
beamdiv=1e-5; %beam divergence approximately
Spotsize=10e-4; %for gaussian dist
G=7;
gammaBeta = sqrt(G^2-1);
dG=1e-6; %energy spread rms control parameter.
zsig=3e-12; %bunch length in seconds
Qbeam=2e6*qe*0;
nps=50000; %number of particles
infile = 'Tomography_test.in';
gdffile = 'results.gdf';
path_to_gpt='C:\Users\vgwr6\OneDrive\Desktop\UCLA\General_Particle_Tracer\bin\';
path_to_in_out_file=' D:\UCLA\PEGASUS\tomography\Generate_dataset\' ; %Must adapt this to your pc

s_beam = [1e-6 0 0 0; 0 2e-8 0 0; 0 0 3e-6 0; 0 0 0 0.5e-8];
centroid_beam = [0 0 0 0];

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
    
    quadvalues = test_settings(1,:)
    I4 = quadvalues(1);
    I5 = quadvalues(2);
    I6 = quadvalues(3);
    fileID = fopen('input_distribution.txt','w');
    fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s\n','x','y','z','Bx','By','Bz','nmacro');
    for i = 1:nps
    fprintf(fileID,'%16.12e %16.12e %16.12e %16.12e %16.12e %16.12e %16.12e\n', x(i), y(i), z(i), Bx(i), By(i), Bz(i), nmacro(i));
    end
    fclose(fileID);
    system(strcat(path_to_gpt,'gpt -o results.gdf ',path_to_in_out_file, 'tomography_test_fixed.in I4=',num2str(I4),' I5=',num2str(I5),' I6=',num2str(I6),' GPTLICENSE=1476385047'),'-echo');
%%
system(strcat(path_to_gpt,'gpt -o results.gdf ', path_to_in_out_file,'tomography_test_fixed.in I4=',num2str(I4),' I5=',num2str(I5),' I6=',num2str(I6),' GPTLICENSE=1476385047'),'-echo');
data=load_gdf(strcat(path_to_in_out_file,gdffile)) ; 
stats=load_gdf(strcat(path_to_in_out_file,'stats_result.gdf')); 
%%
% gammaBeta = mean(phsp(:,6));
dbs = 1.1;
delete_index = [];
for iq = 1:size(test_settings)
    iq
    I4=test_settings(iq,1);
    I5=test_settings(iq,2);
    I6=test_settings(iq,3);
    [I4, I5, I6]
    [R,R_all,z_all,pos] = transport_matrix_ff_focusing_green_var(-[I4, I5, I6],gammaBeta,dbs);
    screen=1; 
    xGPT=data(screen).d.x;
    xpGPT=data(screen).d.Bx;
    yGPT=data(screen).d.y;
    ypGPT=data(screen).d.By;
    phspGPT=zeros(length(x),4);
    phspGPT(:,1)=xGPT;
    phspGPT(:,2)=xpGPT;
    phspGPT(:,3)=yGPT;
    phspGPT(:,4)=ypGPT;
    sigma0=cov(phspGPT);

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
    stdx_macro=zeros(1,length(R_all));
    stdy_macro=zeros(1,length(R_all));

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
    
    if stdx(end)>0.006 || stdy(end)>0.006
        delete_index(end+1) = iq;
    end
end
test_settings(delete_index,:) = [];

%%
x_rot = zeros(size(test_settings(:,1)'));
y_rot = zeros(size(test_settings(:,1)'));
for iq = 1:size(test_settings)
    iq
    I4=test_settings(iq,1);
    I5=test_settings(iq,2);
    I6=test_settings(iq,3);
    quads = -[I4, I5, I6];
    [shearx, sheary, e1, e3, rotationsx, rotationsy] = SER_decomposition(quads, gammaBeta, dbs);
    x_rot(iq) = rotationsx;
    y_rot(iq) = rotationsy;
end
test_settings(:,4) = x_rot;
test_settings(:,5) = y_rot;
%%
figure()
scatter(x_rot, y_rot, 5, 'filled')
xlabel('X Rotation (rad)');
ylabel('Y Rotation (rad)');
%%
sig_x = zeros(size(test_settings(:,1)'));
sig_y = zeros(size(test_settings(:,1)'));
for iq = 1:size(test_settings)
    iq
    %%
    I4=test_settings(iq,1);
    I5=test_settings(iq,2);
    I6=test_settings(iq,3);
    [I4, I5, I6]
    %%
    [R,R_all,z_all,pos] = transport_matrix_ff_focusing_green_var(-[0.06, -0.02, 0],gammaBeta,dbs);
    screen=1; 
    xGPT=data(screen).d.x;
    xpGPT=data(screen).d.Bx;
    yGPT=data(screen).d.y;
    ypGPT=data(screen).d.By;
    phspGPT=zeros(length(x),4);
    phspGPT(:,1)=xGPT;
    phspGPT(:,2)=xpGPT;
    phspGPT(:,3)=yGPT;
    phspGPT(:,4)=ypGPT;
    sigma0=cov(phspGPT);

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
    stdx_macro=zeros(1,length(R_all));
    stdy_macro=zeros(1,length(R_all));

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
    %%
    sig_x(iq) = stdx(end);
    sig_y(iq) = stdy(end);
end
%%
grid_spacing = linspace(0.75, pi, 20);
[grid_x, grid_y] = meshgrid(grid_spacing, grid_spacing);
points = [grid_x(:), grid_y(:)];
%%
rot_pt(:,1) = test_settings(:,4);
rot_pt(:,2) = test_settings(:,5);
[k,dist] = dsearchn(rot_pt, points);
k = unique(k);
rot_grid(:,1) = rot_pt(k,1);
rot_grid(:,2) = rot_pt(k,2);
size(rot_grid)
figure()
scatter(rot_grid(:,1), rot_grid(:,2), 5, 'filled')
xlim([0 pi]);
ylim([0 pi]);
%%
current_fin = test_settings(k,:);
figure()
scatter(current_fin(:,4), current_fin(:,5), 5, 'filled')
xlim([0 pi]);
ylim([0 pi]);
%%
writematrix(current_fin, 'current_settings_G17.csv');
%%
phsp1 = [data(1).d.x data(1).d.Bx data(1).d.y data(1).d.By];
figure()
sgtitle('Initial Distribution')
subplot(2,2,1)
plot(phsp1(:,1),phsp1(:,3),'.')
title('x vs. y')
subplot(2,2,2)
plot(phsp1(:,2),phsp1(:,4),'.')
title('xp vs. yp')
subplot(2,2,3)
plot(phsp1(:,1),phsp1(:,2),'.')
title('x vs. xp')
subplot(2,2,4)
plot(phsp1(:,3),phsp1(:,4),'.')
title('y vs. yp')