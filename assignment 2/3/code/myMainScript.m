seed = 0;
threshold = 1;
dimensions = [256 256];
patch_size = 8;
RMSE_info = "";
images = ["barbara256" "goldhill"];
bases = ["dct2D" "haarWavelet"]
create_directory("../results");
for image = images
    RMSE_info = process(image, bases, dimensions, patch_size, seed, threshold, RMSE_info);
end
file = fopen('../results/RMSE.txt', 'wt');
fprintf(file, RMSE_info);
fclose(file);
function RMSE_info = process(image_name, bases, dimensions, patch_size, seed, threshold, RMSE_info)
    rng(seed);

    image_input = imread("../results/"+ image_name +".png");
    image_input = double(image_input);
    number_of_rows = dimensions(1);
    number_of_columns = dimensions(2);
    image_input = image_input(1:number_of_rows, 1:number_of_columns);
    %crops the loaded image to specified rows and columns
    %image_input variable is reassigned to a submatrix of itself
    image_with_gaussian_noise = image_input + generate_gaussian_noise(dimensions, 0, 4);
    filename = image_name + " with noise";
    save_image(image_with_gaussian_noise, "../results/" + filename + ".png");
    RMSE = calculate_RMSE(image_input, image_with_gaussian_noise);
    RMSE_info = RMSE_info + "RMSE = " + string(RMSE) + " for " + filename + newline;

    measurement_matrix = generate_gaussian_noise([1/2*patch_size*patch_size patch_size*patch_size], 0, 1);
    % 32*64 MeasurementMatrix whose entries are sampled from standard Gaussian
    for basis = bases
        sparse_basis_matrix = feval(basis, patch_size);
        RMSE_info = generate_result("reconstructed using all measurements, without noise",        basis, image_name, image_input,               patch_size, eye(patch_size*patch_size), sparse_basis_matrix, threshold, RMSE_info);
        RMSE_info = generate_result("reconstructed using all measurements, with noise",           basis, image_name, image_with_gaussian_noise, patch_size, eye(patch_size*patch_size), sparse_basis_matrix, threshold, RMSE_info);
        RMSE_info = generate_result("reconstructed using compressive measurements, without noise",basis, image_name, image_input,               patch_size, measurement_matrix,         sparse_basis_matrix, threshold, RMSE_info);
        RMSE_info = generate_result("reconstructed using compressive measurements, with noise",   basis, image_name, image_with_gaussian_noise, patch_size, measurement_matrix,         sparse_basis_matrix, threshold, RMSE_info);
    end
end

function RMSE_info = generate_result(comment, basis,  image_name, image_input, patch_size, measurement_matrix, sparse_basis_matrix, threshold, RMSE_info)
    [image_reconstructed, RMSE] = patch_reconstruct(image_input, patch_size, measurement_matrix, sparse_basis_matrix, threshold);
    create_directory("../results/"+basis);
    save_image(image_reconstructed, "../results/" + basis + "/" + image_name + " " + comment + ".png");
    RMSE_info = RMSE_info + "RMSE = " + string(RMSE) + " for " + image_name + " " + comment + " and " + basis + " basis " + newline;
end

function [image_reconstructed, RMSE] = patch_reconstruct(image_input, patch_size, measurement_matrix, sparse_basis_matrix, threshold)
    dimensions = size(image_input);
    number_of_rows = dimensions(1);
    number_of_columns = dimensions(2);
    number_of_overlapping_patches = zeros(dimensions);                
    % Matrix of size dimensions whose entries are all 0
    image_reconstructed = zeros([number_of_rows, number_of_columns]);
    for i = 1:number_of_rows+1-patch_size
        for j = 1:number_of_columns+1-patch_size
            % disp(string(i)+', '+string(j));
            x = image_input(i:i+patch_size-1, j:j+patch_size-1);
            x = vectorify(x);
            y = measurement_matrix * x;
            theta_hat = ista(y,  measurement_matrix*sparse_basis_matrix, threshold);
            x_hat = sparse_basis_matrix * theta_hat;
            image_reconstructed(i:i+patch_size-1, j:j+patch_size-1) = image_reconstructed(i:i+patch_size-1, j:j+patch_size-1) + reshape(x_hat, [patch_size, patch_size]);
            number_of_overlapping_patches(i:i+patch_size-1, j:j+patch_size-1) = number_of_overlapping_patches(i:i+patch_size-1, j:j+patch_size-1) + ones(patch_size, patch_size);
        end
    end
    image_reconstructed = image_reconstructed ./ number_of_overlapping_patches;
    RMSE = calculate_RMSE(image_input, image_reconstructed);
end

function gaussian_noise_matrix = generate_gaussian_noise(size, mean, variance)
    gaussian_noise_matrix = sqrt(variance)*(mean + randn(size));
    % randn(size):samples 'size' elements from standar gaussian 
    % shifts these to N(mean,variance)
end

function image_vector = vectorify(image)
    image_vector = reshape(image, [prod(size(image)) 1]);
    % convert the image into a matrix of size [prod(size(image)) 1]
end

% Utility functions from my CS663 assignments
function RMSE = calculate_RMSE(image_original, image_reconstructed)
    image_original = vectorify(image_original);
    image_reconstructed = vectorify(image_reconstructed);
    RMSE = norm(image_reconstructed - image_original) / norm(image_original);
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

function [] = create_directory(path)
    if ~exist(path)
        mkdir(path)
    end
end