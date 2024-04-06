seed = 0;
epsilon = 1e-2;
delta_k = 2;
lambdas = [25 50 75 100];
dimensions = [128 128];
f = 0.7;
fraction_validation = 0.9;
image_names = ["textture_sand", "texture_brick"];
size_patch = 7;
size_neighbourhood = 19;

for image_name = image_names
	disp('image = ' + image_name);
	images_reconstructed = cell(1, length(lambdas));
	RMSEs = zeros(1, length(lambdas));
	validation_errors = zeros(1, length(lambdas));
	for i = 1:length(lambdas) % replace "parfor" with "for" if you don't wish to install parallel computing toolbox
		disp(string(lambdas(i)));
		[images_reconstructed{i}, validation_errors(i), RMSEs(i)] = process(seed, image_name, dimensions, f, fraction_validation, lambdas(i), epsilon, delta_k, size_patch, size_neighbourhood);
	end
	[optimal_validation_errors, optimal_lambda_index] = min(validation_errors);
	optimal_RMSE = RMSEs(optimal_lambda_index);
	save("../../media/Q3 "+ image_name + " RMSEs.mat", "RMSEs");
	save_image(images_reconstructed{optimal_lambda_index}, "../../media/Q3 " + image_name + " reconstructed.png")
	disp('RMSE = ' + string(optimal_RMSE) + ' for image ' + image_name + ', with lambda = ' + string(lambdas(optimal_lambda_index)));
end

function [image_reconstructed_measurements, validation_error, RMSE] = process(seed, image_name, dimensions, f, fraction_validation, lambda, epsilon, delta_k, size_patch, size_neighbourhood)
	rng(seed);
	image = imread("../../media/"+ image_name +".tiff");
    image = double(image);
    image = image/255;
    number_of_rows = dimensions(1);
    number_of_columns = dimensions(2);
    image = image(1:number_of_rows, 1:number_of_columns);
	save_image(255*image, "../../media/Q3 " + image_name + " cropped.png")

    m = floor(f*number_of_rows*number_of_columns);
    v = floor(fraction_validation*m);
	mask_measurements = zeros(dimensions);
	mask_validation = zeros(dimensions);
	mask_reconstruction = zeros(dimensions);
	mask = zeros(dimensions);
	indices_measurement = randperm(number_of_rows*number_of_columns, m);
	indices_validation = indices_measurement(1:v);
	indices_reconstruction = indices_measurement(v+1:end);
	mask_measurements(indices_measurement) = 1;
	mask_validation(indices_validation) = 1;
	mask_reconstruction(indices_reconstruction) = 1;

    image_measurements = mask_measurements .* image;
    image_validation = mask_validation .* image;
    save_image(255*image_measurements, "../../media/Q3 " + image_name + " all measurements.png")
    save_image(255*image_validation, "../../media/Q3 " + image_name + " validation.png")
	image_reconstructed_measurements = zeros(size(image));
	image_reconstructed_validation = zeros(size(image));
	number_of_overlapping_patches = zeros(size(image));

	neighbourhood_threshold = floor((size_neighbourhood-1)/2);
	neighbourhood_patch = floor((size_patch-1)/2);
	for i = 1:number_of_rows-size_patch+1
		for j = 1:number_of_columns-size_patch+1
			% disp(string(i)+', '+string(j));
			number_of_overlapping_patches(i:i+size_patch-1, j:j+size_patch-1) = number_of_overlapping_patches(i:i+size_patch-1, j:j+size_patch-1) + 1;
			[left, right, top, bottom] = calculate_neighbourhood_corners(i+neighbourhood_patch, j+neighbourhood_patch, number_of_rows, number_of_columns, neighbourhood_threshold);
			neighbourhood_reconstructed_measurements = process_neighbourhood(image_measurements(top:bottom, left:right), mask_measurements(top:bottom, left:right), size_patch, lambda, epsilon, delta_k);
			neighbourhood_reconstructed_validation = process_neighbourhood(image_measurements(top:bottom, left:right), mask_validation(top:bottom, left:right), size_patch, lambda, epsilon, delta_k);
			
			x = i + neighbourhood_patch - top + 1;
			y = j + neighbourhood_patch - left + 1;
			
			image_reconstructed_measurements(i:i+size_patch-1, j:j+size_patch-1) = image_reconstructed_measurements(i:i+size_patch-1, j:j+size_patch-1) + neighbourhood_reconstructed_measurements(x-neighbourhood_patch:x+neighbourhood_patch, y-neighbourhood_patch:y+neighbourhood_patch);
			image_reconstructed_validation(i:i+size_patch-1, j:j+size_patch-1) = image_reconstructed_validation(i:i+size_patch-1, j:j+size_patch-1) + neighbourhood_reconstructed_validation(x-neighbourhood_patch:x+neighbourhood_patch, y-neighbourhood_patch:y+neighbourhood_patch);
		end
	end
	image_reconstructed_measurements = image_reconstructed_measurements ./ number_of_overlapping_patches;
	RMSE = calculate_RMSE(image, image_reconstructed_measurements);
	image_reconstructed_measurements = image_reconstructed_measurements * 255;

	image_reconstructed_validation = image_reconstructed_validation ./ number_of_overlapping_patches;
	validation_error = mean(mask_reconstruction .* (image - image_reconstructed_validation).^2, "all");
	image_reconstructed_validation = image_reconstructed_validation * 255;
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