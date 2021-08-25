gammaBeta = mean(phsp_clone(:,6));
dbs = 1.1;  
[R,R_all,z_all,pos] = transport_matrix_ff_focusing_green_var(-quadvalues,gammaBeta,dbs);
screen=1; 
x=phsp_clone(:,1);
xp=phsp_clone(:,2);
y=phsp_clone(:,3);
yp=phsp_clone(:,4);
sigma0=cov(phsp_clone(:,1:4))

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
    
figure()
sgtitle('after N iterations')
subplot(2,2,1)
scatter(phsp_clone(:,1),phsp_clone(:,3),10,'filled', 'MarkerFaceAlpha', 0.05)
xlim([-0.01 0.01]);
ylim([-0.01 0.01]);
title('x vs. y')
subplot(2,2,2)
scatter(phsp_clone(:,2),phsp_clone(:,4),10,'filled', 'MarkerFaceAlpha', 0.05)
title('xp vs. yp')
subplot(2,2,3)
scatter(phsp_clone(:,1),phsp_clone(:,2),10,'filled', 'MarkerFaceAlpha', 0.05)
title('x vs. xp')
subplot(2,2,4)
scatter(phsp_clone(:,3),phsp_clone(:,4),10,'filled', 'MarkerFaceAlpha', 0.05)
title('y vs. yp')

   