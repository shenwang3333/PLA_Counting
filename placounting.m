function [res] = placounting (im1, im2)

%% load 2* single-channel .tif images
file_list = dir(strcat('*.tif')); 

im1 = imread(strcat(file_list(1).name));
im2 = imread(strcat(file_list(2).name));


%% cell-segmentation: process nucleus channel [im1] 

% Otsu adaptive binarize, for DAPI/Hoechst-stained nucleus
im1_bw = imbinarize(im1, 'global');
% remove small area, size = s1 (pixel)
% using small size s1 = 1000
s1 = 1000;
im1_bw_t1 = bwareaopen(im1_bw, s1);
% fill holes
im1_bw_t2 = imfill(im1_bw_t1, 'holes');

% construct structral elements, disk diameter = d1
% using disk diameter d1=10
d1 = 10;
SE1 = strel('disk', d1);
% close 
im1_bw_cl = imclose(im1_bw_t2, SE1);

% watershed segmentation, return dist-labeled matrix[seg_Ld]
D_im1 = -bwdist(~im1_bw_cl);
mask = imextendedmin(D_im1, 2);
D_im1_2 = imimposemin(D_im1, mask);

seg_Ld = watershed(D_im1_2);


%% PLA-puncta conting: process PLA channel [im2]
%% use spotdet() function

% bp-filter, ns = 1, os = 3, threshold = 50
im2_bp = bpfilter(im2, 1, 3, 50);

% detect spots
% using threshold = 5 (as bpfilter transform the data to double); os_size = 7
thresh = 2;
os_size = 7;

im2_det = spotdet(im2_bp, thresh, os_size)

% projection of spot center to spot matrix
[row, col] = size(im2);
spot_mat = zeros(row, col);

spot_mat(sub2ind(size(spot_mat), im2_det(:, 1), im2_det(:, 2))) = 1;

% convert to binary
spot_mat_bw = imbinarize(spot_mat);


%% counting PLA puncta per cell

% counting region number
cell_num = max(max(seg_Ld));

% results matrix
res = zeros(cell_num, 1);

for j = 1:cell_num

	cell_mat = zeros(row, col);

	[row2, col2] = find(seg_Ld == j);

	cell_mat(sub2ind(size(cell_mat), row2, col2)) = 1;

	cell_mat_bw = imbinarize(cell_mat);

	spot_num_mat = spot_mat_bw & cell_mat_bw;

	PLA_counts = sum(sum(spot_num_mat));

	res(j, 1) = PLA_counts

end




