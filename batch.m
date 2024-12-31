
clear;

folder_list = dir(strcat('*.frames'));
folder_num = length(folder_list);

final_res = {};

for j = 1:folder_num

	cd(strcat(folder_list(j).name));

	res_temp = placounting;

	final_res{j, 1} = res_temp;

	cd ..

end

final_res_array = cell2mat(final_res);

