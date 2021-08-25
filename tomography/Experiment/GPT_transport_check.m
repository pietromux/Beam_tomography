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
nps=50000; %number of particles

s_beam = [1e-6 0 0 0; 0 2e-8 0 0; 0 0 3e-6 0; 0 0 0 0.5e-8];
centroid_beam = [0 0 0 0];
r0p = initial_momentum_guess();
generate_initial_phsp;
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
    path_to_gpt='C:\Users\vgwr6\OneDrive\Desktop\UCLA\General_Particle_Tracer\bin\' ;
    fileID = fopen('input_distribution.txt','w');
    fprintf(fileID,'%6s %6s %6s %6s %6s %6s %6s\n','x','y','z','Bx','By','Bz','nmacro');
    for i = 1:nps
    fprintf(fileID,'%16.12e %16.12e %16.12e %16.12e %16.12e %16.12e %16.12e\n', x(i), y(i), z(i), Bx(i), By(i), Bz(i), nmacro(i));
    end
    fclose(fileID);
    system(strcat(path_to_gpt,'asci2gdf -o start.gdf input_distribution.txt'));
    
    %%
    quadvalues = quadlist(38,:)
    
    I4=quadvalues(1) ;
    I5=quadvalues(2);
    I6=quadvalues(3);
    %%
%     quadvalues=[I4 I5 I6]
%     quadvalues = [-0.0600    0.0200         0];
    infile = 'Tomography_test.in';
    gdffile = 'results.gdf';
    path_to_gpt='C:\Users\vgwr6\OneDrive\Desktop\UCLA\General_Particle_Tracer\bin\';
    path_to_in_out_file=' D:\UCLA\PEGASUS\tomography\Generate_dataset\' ; %Must adapt this to your pc

    system(strcat(path_to_gpt,'gpt -v -o',path_to_in_out_file,'results.gdf ',path_to_in_out_file, 'Tomography_test.in',' I4=',num2str(I4),' I5=',num2str(I5),' I6=',num2str(I6),' GPTLICENSE=1476385047'),'-echo');
    %%
    data=load_gdf(strcat(path_to_in_out_file,gdffile)) ; 
    callGPT2=strcat(path_to_gpt,'gdfa -o',path_to_in_out_file,'stats_result.gdf ',path_to_in_out_file,gdffile,' position ',' stdBy ',' stdBx ',' stdBz',' stdx ',' stdy',' CSgammax',' CSgammay',' CSalphax',' CSalphay', ' CSbetax',' CSbetay',' nemixrms',' nemiyrms',' nemirrms',' nemizrms',' avgz',' stdz',' avgG',' stdG');
    system(callGPT2,'-echo') ;
    stats=load_gdf(strcat(path_to_in_out_file,'stats_result.gdf')); 
    %% compare envlope evolution
    
    gammaBeta = mean(phsp(:,6));
    dbs = 1.1;  
    [R,R_all,z_all,pos] = transport_matrix_ff_focusing_green_var(-quadvalues,gammaBeta,dbs);
    screen=1; 
    xGPT=data(screen).d.x;
    xpGPT=data(screen).d.Bx;
    yGPT=data(screen).d.y;
    ypGPT=data(screen).d.By;
    phspGPT=zeros(length(xGPT),4);
    phspGPT(:,1)=xGPT;
    phspGPT(:,2)=xpGPT;
    phspGPT(:,3)=yGPT;
    phspGPT(:,4)=ypGPT;
    sigma0=cov(phspGPT)

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
    
    

%% compare GPT to the matrix code
figure()
plot(z_all,stdx,'LineWidth',2)
hold on
plot(z_all,stdy,'LineWidth',2)
hold on
plot(stats.d.position(2:end),stats.d.stdx(2:end),'--','LineWidth',2)
hold on
plot(stats.d.position(2:end),stats.d.stdy(2:end),'--','LineWidth',2)
legend('\sigma_x(z)','\sigma_y(z)','\sigma_x(z) (GPT)','\sigma_y(z) (GPT)','Location','best')
xlabel('z(m)')
set(gca,'FontSize',15)

%%
figure()
sgtitle('actual initial phase space')
subplot(2,2,1)
plot(xGPT,yGPT,'.')
title('x vs. y')
subplot(2,2,2)
plot(xpGPT,ypGPT,'.')
title('xp vs. yp')
subplot(2,2,3)
plot(xGPT,xpGPT,'.')
title('x vs. xp')
subplot(2,2,4)
plot(yGPT,ypGPT,'.')
title('y vs. yp')
%%
figure()
sgtitle('actual initial phase space')
subplot(2,2,1)
dscatter(xGPT,yGPT)
xlim([1.1*min(min(xGPT)) 1.1*max(max(xGPT))]);
ylim([1.1*min(min(yGPT)) 1.1*max(max(yGPT))]);
colorbar
title('x vs. y')
subplot(2,2,2)
dscatter(xpGPT,ypGPT)
xlim([1.1*min(min(xpGPT)) 1.1*max(max(xpGPT))]);
ylim([1.1*min(min(ypGPT)) 1.1*max(max(ypGPT))]);
colorbar
title('xp vs. yp')
subplot(2,2,3)
dscatter(xGPT,xpGPT)
xlim([1.1*min(min(xGPT)) 1.1*max(max(xGPT))]);
ylim([1.1*min(min(xpGPT)) 1.1*max(max(xpGPT))]);
colorbar
title('x vs. xp')
subplot(2,2,4)
dscatter(yGPT,ypGPT)
xlim([1.1*min(min(yGPT)) 1.1*max(max(yGPT))]);
ylim([1.1*min(min(ypGPT)) 1.1*max(max(ypGPT))]);
colorbar
title('y vs. yp')