clear all
close all

%% parse input filenames and extract quadrupole values
inputdir = strcat(pwd,'\Generate_dataset\dg_gpt\');
D = dir(inputdir);
p = 1;
for jf = 1:size(D,1)
    file_image(jf) = contains(D(jf).name,'.bmp') & contains(D(jf).name,'target') ;
    if (file_image(jf))
    filenamelist(p) = string(D(jf).name);
    quadv = split(filenamelist(p),'_');
    I4 = str2num(quadv(2));
    I5 = str2num(quadv(3));
    auxI6 = split(quadv(4),'.bmp');
    I6 = str2num(auxI6(1));
    quadlist(p,:) = [I4, I5, I6];
    p = p+1;
    end
end

%% Initial phase space guess
%r0p = initial_momentum_guess();
r0p = 7e-3;

%%
generate_initial_phsp;
npart = size(phsp,1);


%% Iterate over quadrupole values
tic
iteration = 0;

%loop_size = 1;
%conv = npart * size(filenamelist,2);
%reassign_count = size(filenamelist,2) * 20;
%respawn_idx = zeros(1, npart) + reassign_count ;

phsp0 = phsp;
nrespawn = npart;
idx = ones(1,npart);

while nrespawn > 0.05*npart && iteration <= 150
    iteration = iteration + 1
%    reassign_count = 0;
%    reassign_ratio = 0;
%    phsp_tot = zeros(size(phsp(:,1:2)))
%    conv = 0;

    tot_idx = zeros(1, npart);
    tot_weights = zeros(1,npart);
    
    % add comparison with scr4 figure
    target = B/sum(sum(B))*npart;
    phsp_f = phsp;
    [fix_idx, neg_x, neg_y, weights, conv_param] = compare_images(phsp, phsp_f, target);
    tot_idx = tot_idx + 5*fix_idx;
    tot_weights = tot_weights + 5*weights;
    
    for nq = 1:size(filenamelist,2);
        quadvalues = quadlist(nq,:);
        filename = filenamelist(nq);
        target = double(imread(strcat(inputdir,filename)));  
        norm = sum(sum(target));
        target = target/norm*npart;       
        
        %% transport and compare images
        [phsp_f, R] = matrix_quad_transport(phsp,quadvalues);
        [fix_idx, neg_x, neg_y, weights, conv_param] = compare_images(phsp, phsp_f, target);
        tot_idx = tot_idx + fix_idx;
        tot_weights = tot_weights +weights;
    end
%         conv_hist = [conv_param];
%         while abs(conv_hist(end)) > 0.05 * npart
%             [phsp_f, R] = matrix_quad_transport(phsp_clone,quadvalues);
%             [idx, neg_x, neg_y, weights, ~] = compare_images(phsp_clone, phsp_f, target);
%             nrespawn_init = sum(idx);
%             phsp_re = momentum_reassign(phsp_clone, neg_x, neg_y, weights, idx,r0p,'gaussian', R);
%             [phsp_f, R] = matrix_quad_transport(phsp_re,quadvalues);
%             [idx, neg_x, neg_y, weights, conv_param] = compare_images(phsp_re, phsp_f, target);
%             nrespawn = sum(idx)
%             conv_part = nrespawn_init - nrespawn;
%             conv_hist(end+1) = conv_param
%             if length(conv_hist)>5
%                 break
%             elseif abs(conv_hist(end)) > abs(conv_hist(end-1))
%                 break
%             else
%                 phsp_clone = phsp_re;
%             end
%         end
%         conv_ratio = sum(nrespawn_init, 'all') / npart;
%         reassign_ratio = reassign_ratio + conv_ratio;
%             msp(:,1) = phsp_clone(:,2);
%             msp(:,2) = phsp_clone(:,4);
%             phsp_tot = phsp_tot + msp * conv_ratio;
%         reassign_count = reassign_count + 1
%         respawn_idx = respawn_idx + fix_idx;
%         conv = conv + conv_hist(end);   
%    end
% kept_idx = (respawn_idx < (0.05 * reassign_count));
% sum(kept_idx, 'all');
% sum(respawn_idx, 'all')
% phsp(:,2) = phsp(:,2).* kept_idx' + ~kept_idx' .* phsp_tot(:,1)/ reassign_ratio;
% phsp(:,4) = phsp(:,4).* kept_idx' + ~kept_idx' .* phsp_tot(:,2)/ reassign_ratio;

%% 3- reassign momentum - choose max number of images to miss. 10 in this example
        idx = idx & (tot_idx> 2);
        figure(10)
        plot(tot_weights, tot_idx,'.')
        nrespawn = sum(idx);
        [nrespawn, npart-nrespawn]
        phsp_re = momentum_reassign_simple(phsp,idx,r0p,'fit');
        phsp = phsp_re;
toc
%% visualize phase space at N iteration
figure(5)
titl = strcat('after',{' '}, num2str(iteration),' iterations');
sgtitle(titl)
subplot(2,2,1)
dscatter(phsp(:,1),phsp(:,3))
xlim([1.1*min(min(phsp0(:,1))) 1.1*max(max(phsp0(:,1)))]);
ylim([1.1*min(min(phsp0(:,3))) 1.1*max(max(phsp0(:,3)))]);
title('x vs. y')
colorbar
subplot(2,2,2)
dscatter(phsp(:,2),phsp(:,4))
xlim([1.1*min(min(phsp0(:,2))) 1.1*max(max(phsp0(:,2)))]);
ylim([1.1*min(min(phsp0(:,4))) 1.1*max(max(phsp0(:,4)))]);
title('xp vs. yp')
colorbar
subplot(2,2,3)
dscatter(phsp(:,1),phsp(:,2))
xlim([1.1*min(min(phsp0(:,1))) 1.1*max(max(phsp0(:,1)))]);
ylim([1.1*min(min(phsp0(:,2))) 1.1*max(max(phsp0(:,2)))]);
title('x vs. xp')
colorbar
subplot(2,2,4)
dscatter(phsp(:,3),phsp(:,4))
xlim([1.1*min(min(phsp0(:,3))) 1.1*max(max(phsp0(:,3)))]);
ylim([1.1*min(min(phsp0(:,4))) 1.1*max(max(phsp0(:,4)))]);
title('y vs. yp')
colorbar
end
s_beam_rec = cov([phsp(:,1) phsp(:,2) phsp(:,3) phsp(:,4)]);
s_beam_target
s_beam_rec