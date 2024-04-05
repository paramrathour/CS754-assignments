seed = 0;
epsilon = 0.01;
delta_k = 2;
% lambdas = 10 .^ (1:3);
lambda = 100;
dimensions = [128 128];
dimensions = [64 64];
f = 0.7;
image_names = ["textture_sand", "texture_brick"];
size_patch = 7;
size_neighbourhood = 19;

for image_name = image_names
	RMSE = process(image_name, dimensions, f, seed, lambda, epsilon, delta_k, size_patch, size_neighbourhood)
end

function RMSE = process(image_name, dimensions, f, seed, lambda, epsilon, delta_k, size_patch, size_neighbourhood)
	rng(seed);
	image = imread("../../media/"+ image_name +".tiff");
    image = double(image);
    image = image/255;
    number_of_rows = dimensions(1);
    number_of_columns = dimensions(2);
    image = image(1:number_of_rows, 1:number_of_columns);
	save_image(255*image, "../../media/Q3 " + image_name + " cropped.png")

    m = floor(f*number_of_rows*number_of_columns);
	mask = zeros([number_of_rows number_of_columns]);
	mask(randperm(number_of_rows*number_of_columns, m)) = 1;
    image_measured = mask .* image;
    save_image(255*image_measured, "../../media/Q3 " + image_name + " measured.png")
	image_reconstructed = zeros(size(image));
	number_of_overlapping_patches = zeros(size(image));

	neighbourhood_threshold = floor((size_neighbourhood-1)/2);
	neighbourhood_patch = floor((size_patch-1)/2);
	for i = 1+neighbourhood_patch:number_of_rows-neighbourhood_patch
		for j = 1+neighbourhood_patch:number_of_columns-neighbourhood_patch
			disp(string(i)+', '+string(j));
			number_of_overlapping_patches(i-neighbourhood_patch:i+neighbourhood_patch, j-neighbourhood_patch:j+neighbourhood_patch) = number_of_overlapping_patches(i-neighbourhood_patch:i+neighbourhood_patch, j-neighbourhood_patch:j+neighbourhood_patch) + 1;
			[left, right, top, bottom] = calculate_neighbourhood_corners(i, j, number_of_rows, number_of_columns, neighbourhood_threshold);
			neighbourhood_reconstructed = process_neighbourhood(image_measured(left:right, top:bottom), mask(left:right, top:bottom), size_patch, lambda, epsilon, delta_k);
			x = i - top + 1;
			y = j - left + 1;
			image_reconstructed(i-neighbourhood_patch:i+neighbourhood_patch, j-neighbourhood_patch:j+neighbourhood_patch) = image_reconstructed(i-neighbourhood_patch:i+neighbourhood_patch, j-neighbourhood_patch:j+neighbourhood_patch) + neighbourhood_reconstructed(x-neighbourhood_patch:x+neighbourhood_patch, y-neighbourhood_patch:y+neighbourhood_patch);
		end
	end
	image_reconstructed = image_reconstructed ./ number_of_overlapping_patches;
	RMSE = calculate_RMSE(image, image_reconstructed);
	save_image(255*image_reconstructed, "../../media/Q3 " + image_name + " reconstructed.png")
end

function neighbourhood_reconstructed = process_neighbourhood(neighbourhood, mask, size_patch, lambda, epsilon, delta_k)
	number_of_rows = size(neighbourhood, 1);
	number_of_columns = size(neighbourhood, 2);
	patch_matrix = zeros([size_patch*size_patch (number_of_rows+1-size_patch)*(number_of_columns+1-size_patch)]);
	mask_matrix = zeros([size_patch*size_patch (number_of_rows+1-size_patch)*(number_of_columns+1-size_patch)]);
	counter = 1;
	[X, Y] = ndgrid(1:number_of_rows-size_patch+1, 1:number_of_columns-size_patch+1);
	X = X(:);
	Y = Y(:);
	neighbourhood_reconstructed = zeros(size(neighbourhood));
	number_of_overlapping_patches = zeros(size(neighbourhood));
	for k = 1:length(X)
		x = X(k); y = Y(k);
		number_of_overlapping_patches(x:x+size_patch-1, y:y+size_patch-1) = number_of_overlapping_patches(x:x+size_patch-1, y:y+size_patch-1) + 1;
		patch_image = neighbourhood(x:x+size_patch-1, y:y+size_patch-1);
		patch_mask = mask(x:x+size_patch-1, y:y+size_patch-1);
		patch_matrix(:, k) = patch_image(:);
		mask_matrix(:, k) = patch_mask(:);
	end
	patches_reconstructed = SVT(patch_matrix, mask_matrix, lambda, epsilon, delta_k);
	for k = 1:length(X)
		x = X(k); y = Y(k);
		neighbourhood_reconstructed(x:x+size_patch-1, y:y+size_patch-1) = neighbourhood_reconstructed(x:x+size_patch-1, y:y+size_patch-1) + reshape(patches_reconstructed(:, k), [size_patch, size_patch]);
	end
	neighbourhood_reconstructed = neighbourhood_reconstructed ./ number_of_overlapping_patches;
end

function RMSE = calculate_RMSE(image_original, image_reconstructed)
    RMSE = norm(image_original - image_reconstructed, "fro") / norm(image_original, "fro");
end

% Utility functions from my CS663 assignments
function [left, right, top, bottom] = calculate_neighbourhood_corners(i, j, number_of_rows, number_of_columns, threshold)
    left = max(1, j - threshold);
    right = min(number_of_columns, j + threshold);
    top = max(1, i - threshold);
    bottom = min(number_of_rows, i + threshold);
end

function image_input = image_bound(image_input)
    image_input(image_input < 0) = 0;
    image_input(image_input > 255)= 255;
end

function image_integer = image_double_to_integer(image_double)
    image_integer = uint8(image_double);
end

function [] = save_image(image, filename)
    image = image_bound(image);
    image = image_double_to_integer(image);
    imwrite(image, filename);
end