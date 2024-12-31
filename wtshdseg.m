function [label_wtsd] = wtshdseg(input_img, object_scale)

%%
%%

sz_img = size(input_img);

if numel(sz_img) > 2
	
	% forced convert to monoimage
	img = rgb2gray(input_img);

else

	img = input_img;

end

% calculate raw gradient magnitude of original image
gradmag_raw = imgradient(img);
% calculate raw label matrix with watershed algorithm
label_raw = watershed(gradmag_raw);

% construct structral elements
% for open/close reconstruction
SE1 = strel('disk', floor(object_scale./2));


% open-reconstruction
img_erode = imerode(img, SE1);
img_er_recon = imreconstruct(img_erode, img);

% successive close-reconstruction
img_er_recon_dil = imdilate(img_er_recon, SE1);
img_oc_recon = imcomplement(imreconstruct(imcomplement(img_er_recon_dil),imcomplement(img_er_recon)));


% calculate regional maximal patches
rmax_bw = imregionalmax(img_oc_recon);

% if the objects are too small, it will be processed with minimal structral elements by image open/close operation after image reconstruction 
if floor(object_scale./4) > 2

	SE2 = strel(ones(floor(object_scale./4)));

	rmax_bw2 = imopen(rmax_bw, SE2);

else
	
	SE3 = strel(ones(2));

	rmax_bw2 = imopen(rmax_bw, SE3);

end

% remove small patches
rmax_bw3 = bwareaopen(rmax_bw2, object_scale.*2);

% calculate binary image of the original image processed by open-close reconstruction
% use Otsu algorithm
img_oc_recon_bw = imbinarize(img_oc_recon, 'global');

% calculate distance matrix of the binary image and perform watershed segmentation
dist_mat = bwdist(img_oc_recon_bw);
label_dist_mat = watershed(dist_mat);
% watershed ridge lines generated from image processed by open-close reconstruction
ridge_line_mat = label_dist_mat == 0;

% final watershed segmentation
% modify reginal-min using imimposemin function
gradmag = imimposemin(gradmag_raw, ridge_line_mat | rmax_bw3);

label_wtsd = watershed(gradmag);



