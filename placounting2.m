function [res] = placounting2()

%% load 2* single-channel .tif images
file_list = dir(strcat('*.tif')); 

im1 = imread(strcat(file_list(2).name));
im2 = imread(strcat(file_list(3).name))


%% cell-segmentation: process whole cell channel [im1] 

% convert to binary
% using scale bw_param = 0.01
bw_param = 0.01
im1_bw = im2bw(im1, bw_param);

% open operation
% using SE = 3, disk
SE_bw_size = 3;
SE_bw = strel('disk', SE_bw_size);
im1_bw_op = imopen(im1_bw, SE_bw);

% remove small area, size = s1 (pixel)
% using small size s1 = 250
s1 = 250;
im1_bw_t1 = bwareaopen(im1_bw_op, s1);
% fill holes
im1_bw_t2 = imfill(im1_bw_t1, 'holes');

% convert to numerical image
im1_bw_t2_d = im2double(im1_bw_t2);

% refined watershed segmentation, recall function:[wtshdseg]
% using object scale wtshd_obj_scale = 20
wtshd_obj_scale = 20;
seg_Ld = wtshdseg(im1_bw_t2_d, wtshd_obj_scale);


%% PLA-puncta conting: process PLA channel [im2]
%% use spotdet() function

% bp-filter, ns = 1, os = 3, threshold = NaN
im2_bp = bpfilter(im2, 1, 3);

% detect spots
% using threshold = 2 (as bpfilter transform the data to double); os_size = 7
thresh = 2;
os_size = 7;

im2_det = spotdet(im2_bp, thresh, os_size)

% projection of spot center to spot matrix
[row, col] = size(im2);
total_pixel = row*col;
spot_mat = zeros(row, col);

spot_mat(sub2ind(size(spot_mat), im2_det(:, 1), im2_det(:, 2))) = 1;

% convert to binary
spot_mat_bw = imbinarize(spot_mat);


%% counting PLA puncta per cell

% counting region number
cell_num = max(max(seg_Ld));

% results matrix
res = zeros(cell_num-1, 1);

for j = 1:cell_num

	cell_mat = zeros(row, col);

	[row2, col2] = find(seg_Ld == j);

	% discard background area
	% affirm <1/10 total image pixels as background area 
	[sz_area, ~] = size(row2);

	if sz_area < 0.05.*total_pixel

		cell_mat(sub2ind(size(cell_mat), row2, col2)) = 1;

		cell_mat_bw = imbinarize(cell_mat);

		spot_num_mat = spot_mat_bw & cell_mat_bw;

		PLA_counts = sum(sum(spot_num_mat));

		res(j, 1) = PLA_counts;

	else
		continue;
	end

end




