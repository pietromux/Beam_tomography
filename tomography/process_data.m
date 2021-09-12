filepath = 'D:\UCLA\PEGASUS\dataset_8_19_21\';
data_mat = readmatrix(strcat(filepath, 'tempdat.csv'));
[~, idx] = unique(data_mat(:,1:3), 'rows');
data_u = data_mat(idx,:);
size_im = round([max(max(data_u(:,[7,8]))) max(max(data_u(:,[7,8])))] .* 6)
%%
inputdir = 'D:\UCLA\PEGASUS\dataset_8_19_21\';
D = dir(inputdir);
p = 1;
count = 0;
for id = 1:size(data_u,1)
    I4 = data_u(id,1);
    I5 = data_u(id,2);
    I6 = data_u(id,3);
    fileid = data_u(id,4);
    xc = data_u(id,5);
    yc = data_u(id,6);
    stx = data_u(id,7);
    sty = data_u(id,8);
    xmin = xc - size_im(1)/2;
    ymin = yc - size_im(2)/2;
    for jf = 1:size(D,1)
        file_image(jf) = contains(D(jf).name,'.tiff') & contains(D(jf).name, num2str(fileid)) ;
        if (file_image(jf))
            filenamelist(p) = string(D(jf).name);
            quadlist(p,:) = data_u(id,1:3);
            im = double(imread(strcat(inputdir,filenamelist(p))));
            figure(500)
            imagesc(im)
            colorbar
            im = double(im);
            
            x = [50 100 150  160 160 170 180 200 250 270 275 270 260 250 240 230  190 150 100 50  50];
            y = [80 55  55   55  130  130 55  55  100 135 170 200 240 260 270 280  300 305 300 290 80];
            mask = poly2mask(x,y,size(im,1),size(im,2));
            idx = find(~mask .* im>0);
            im(idx) = 0;
%             figure(501)
%             imagesc(im)
%             colorbar
            crop_im = imcrop(im,[xmin ymin size_im(1) size_im(2)]);
%             figure(101)
%             imagesc(crop_im)
%             colorbar
%             crop_im = uint8(crop_im/sum(crop_im,'all')*50000);
            crop_im = double(crop_im)-1*mean(mean(crop_im));%convert int to double.
            crop_im(crop_im<0)=0;
            crop_im = uint8(crop_im/max(max(crop_im))*255-1);
            
            
%             crop_im = imgaussfilt(crop_im, 1);
            crop_im = imresize(crop_im, [800 800]);
            
%             figure(202)
%             imagesc(crop_im)
%             colorbar
%             if sum(crop_im,'all')>100000
%                 count = count+1
%                 image_name = fullfile(strcat('D:\UCLA\PEGASUS\8_19_21_crop\good_im\', 'target_' ,num2str(I4),'_',num2str(I5),'_',num2str(I6),'.bmp'));
%                 imwrite(crop_im, image_name);
%                 crop_im = imread(image_name);
%                 figure(102)
%                 imagesc(crop_im)
%                 colorbar
%                 p = p+1;
%             end
            image_name = fullfile(strcat('D:\UCLA\PEGASUS\8_19_21_crop\', 'target_' ,num2str(I4),'_',num2str(I5),'_',num2str(I6),'.bmp'));
            imwrite(crop_im, image_name);
            crop_im = imread(image_name);
            figure(102)
            imagesc(crop_im)
            colorbar
            p = p+1;
        end
    end
end
% for id = 1:size(data_u,1)
%     I4 = data_u(id,1);
%     I5 = data_u(id,2);
%     I6 = data_u(id,3);
%     fileid = data_u(id,4);
%     xc = data_u(id,5);
%     yc = data_u(id,6);
%     stx = data_u(id,7);
%     sty = data_u(id,8);
%     xmin = xc - size_im(1)/2;
%     ymin = yc - size_im(2)/2;
%     for jf = 1:size(D,1)
%         file_image(jf) = contains(D(jf).name,'.tiff') & contains(D(jf).name, num2str(fileid)) ;
%         if (file_image(jf))
%             filenamelist(p) = string(D(jf).name);
%             quadlist(p,:) = data_u(id,1:3);
%             im = double(imread(strcat(inputdir,filenamelist(p))));
%             figure(200)
%             imagesc(im)
%             colorbar
% %             im = double(im)-mean(im(im>0));%convert int to double.
% %             im(im<0)=0;
%             im(idx) = 0;
%             figure(201)
%             imagesc(im)
%             colorbar
%             crop_im = imcrop(im,[xmin ymin size_im(1) size_im(2)]);
%             figure(101)
%             imagesc(crop_im)
%             colorbar
%             crop_im = uint8(crop_im/sum(crop_im,'all')*50000);
%             if max(max(crop_im))>255
%                 crop_im = crop_im/max(max(crop_im))*255-1;
%             end
%             crop_im = imgaussfilt(crop_im, 0.5);
%             sum(crop_im,'all')
%             
%             figure(202)
%             imagesc(crop_im)
%             colorbar
%             image_name = fullfile(strcat('D:\UCLA\PEGASUS\8_19_21_crop\', 'target_' ,num2str(I4),'_',num2str(I5),'_',num2str(I6),'.bmp'));
%             imwrite(crop_im, image_name);
%             crop_im = imread(image_name);
%             figure(102)
%             imagesc(crop_im)
%             colorbar
%             p = p+1;
%         end
%     end
% end

