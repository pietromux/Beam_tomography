function [idx, neg_x, neg_y, weights, conv_param] = compare_images(phsp_i, phsp_f, target)

% This function compare the distribution obtained from the guess with the measured image and identifies the particles to be killed and respawn
display = 0;
imsize = size(target);
pixelcal = 40e-6;
psf = 40e-6;

A = create_image_from_partcoord(phsp_f, imsize, pixelcal, psf, 0, '');
A = double(A);
A = A/sum(sum(A))*sum(sum(target));
difference = A-target;

if(display)
figure(101)
subplot(1,2,1)
imagesc(A)
colorbar
c2 = caxis;
title('guess')

subplot(1,2,2)
imagesc(target)
colorbar
c1 = caxis;
c3 = [min([c1 c2]), max([c1 c2])];
caxis(c3)
title('target')

subplot(1,2,1)
caxis(c3)
end

weights = zeros(1,size(phsp_f,1));

% Determine if particle j is in a region of excess or lack of charge
for j = 1:size(phsp_f,1)
    ix = round(phsp_f(j,1)/pixelcal)+imsize(1)/2;
    iy = round(-phsp_f(j,3)/pixelcal)+imsize(2)/2;
    if (ix<1 | ix>imsize(1) | iy<1 | iy>imsize(2))
       idx(j) = true ;
       weights(j) = 1;
    else
       idx(j) = (difference(iy,ix) > 0.8*sqrt(target(iy,ix)));
       % 0.8 depends on the filtering and therefore on the ratio between psf and pixelcal
       % for a pure Poisson statistics, it should be 1
       weights(j) = difference(iy,ix);
    end    
end

% ??
%weights = [];
%linear_idx = find(difference < 0);
%weights = round(abs(difference(linear_idx)));
%[row,col] = ind2sub(imsize,linear_idx);
%neg_x = (col - imsize(2)/2) * pixelcal - rand(1,1)* pixelcal;
%neg_y = -(row - imsize(2)/2) * pixelcal - rand(1,1)* pixelcal;
%neg_x = repelem(neg_x, weights);
%neg_y = repelem(neg_y, weights);
neg_x = 0;
neg_y = 0;

if(display)
figure(102)
subplot(1,2,1)
imagesc(difference)
title('difference')
colorbar

subplot(1,2,2)
plot(phsp_i(idx(:),1),phsp_i(idx(:),3),'.')
title('particles to respawn')
end

% if (display)
% figure(103)
% sgtitle('after N iterations')
% subplot(2,2,1)
% scatter(phsp_i(:,1),phsp_i(:,3),10,'filled', 'MarkerFaceAlpha', 0.05)
% xlim([-0.01 0.01]);
% ylim([-0.01 0.01]);
% title('x vs. y')
% subplot(2,2,2)
% scatter(phsp_i(:,2),phsp_i(:,4),10,'filled', 'MarkerFaceAlpha', 0.05)
% title('xp vs. yp')
% subplot(2,2,3)
% scatter(phsp_i(:,1),phsp_i(:,2),10,'filled', 'MarkerFaceAlpha', 0.05)
% title('x vs. xp')
% subplot(2,2,4)
% scatter(phsp_i(:,3),phsp_i(:,4),10,'filled', 'MarkerFaceAlpha', 0.05)
% title('y vs. yp')
% end

conv_param = sum(difference(difference<0));
end
