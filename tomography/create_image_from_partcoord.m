function imphsp = create_image_from_partcoord(phsp, imsize, pixelcal, psf, saveflag, savefilename)

% This function creates an imsize-large image from a particle coordinate
% list phsp (x, px, y, py). pixelcal is the screen calibration and psf is
% the standard deviation of the gaussian spoint spread function of the screen 

% saveflag is an option to save output image to savefilename

n = imsize(1);

xgrid = ((1:imsize(1))-imsize(1)/2)*pixelcal;
ygrid = ((1:imsize(2))-imsize(2)/2)*pixelcal;
xf = phsp(:,1)';
yf = phsp(:,3)';

xr = interp1(xgrid,1:imsize(1),xf,'nearest')';
yr = interp1(ygrid,1:imsize(2),-yf,'nearest')';
iout = isnan(xr) | isnan(yr);
xr(iout)=[];
yr(iout)=[];

Z = accumarray([yr xr],1,[n n]);
imphsp = imgaussfilt(Z,psf/pixelcal);

% need to be careful before saving and check that the saved image has the correct scale
imphsp_aux = imphsp/max(max(imphsp))*255;
imphsp = uint8(imphsp_aux-1);

if saveflag
     imwrite(imphsp,savefilename);
end

end

