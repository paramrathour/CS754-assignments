dimensions = [256 256];
seed = 0;
patch_size = 8;
images = ["barbara256" "goldhill"];
algorithms = ["ista"];
% bases = ["dct2D" "haarWavelet"];
bases = ["dct2D"]; % don't know how to do haarWavelet yet
create_directory("../results");
for image = images
    for basis = bases
        for algorithm = algorithms
            process(image, algorithm, basis, dimensions, patch_size, seed);
        end
    end
end

function [] = process(image_name, algorithm, basis, dimensions, patch_size, seed)
    rng(seed);

    image_input = imread("../results/"+ image_name +".png");
    image_input = double(image_input);
    number_of_rows = dimensions(1);
    number_of_columns = dimensions(2);
    image_input = image_input(1:number_of_rows, 1:number_of_columns);

    image_with_gaussian_noise = image_input + image_bound(generate_gaussian_noise(dimensions, 0, 4));
    filename = "../results/" + image_name + " with noise" + ".png";
    save_image(image_with_gaussian_noise, filename);

    image_reconstructed = zeros([number_of_rows, number_of_columns]);
    measurement_matrix = generate_gaussian_noise([1/2*patch_size*patch_size patch_size*patch_size], 0, 1);
    sparsifying_matrix = feval(basis, patch_size);
    number_of_overlapping_patches = zeros(size(image_input));
    for i = 1:number_of_rows+1-patch_size
        for j = 1:number_of_columns+1-patch_size
            x = image_with_gaussian_noise(i:i+patch_size-1, j:j+patch_size-1);
            x = reshape(x, [patch_size*patch_size 1]);
            y = measurement_matrix * x;
            theta_hat = feval(algorithm, y,  measurement_matrix*sparsifying_matrix);
            x_hat = sparsifying_matrix * theta_hat;
            image_reconstructed(i:i+patch_size-1, j:j+patch_size-1) = reshape(x_hat, [patch_size, patch_size]);
            number_of_overlapping_patches(i:i+patch_size-1, j:j+patch_size-1) = number_of_overlapping_patches(i:i+patch_size-1, j:j+patch_size-1) + ones(patch_size, patch_size);
        end
    end
    image_reconstructed = image_reconstructed ./ number_of_overlapping_patches;
    filename = "../results/" + image_name + ", basis = " + basis + ", algorithm = " + algorithm + ".png";
    save_image(image_reconstructed, filename);
    
    RMSE = calculate_RMSE(image_input, image_reconstructed, dimensions);
    disp(["RMSE for " + image_name + ", with " + basis + " and " + algorithm + " = " + string(RMSE)]);
end

function gaussian_noise_matrix = generate_gaussian_noise(size, mean, sigma)
    gaussian_noise_matrix = sigma*(mean + randn(size));
end

% Utility functions from my CS663 assignments
function RMSE = calculate_RMSE(image_original, image_reconstructed, dimensions)
    image_original = reshape(image_original, [dimensions(1)*dimensions(2) 1]);
    image_reconstructed = reshape(image_reconstructed, [dimensions(1)*dimensions(2) 1]);
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