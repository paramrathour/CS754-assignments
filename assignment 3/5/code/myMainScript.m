seed = 0;
threshold = 1;
% dimensions = [120 240];
dimensions = [288 352];
patch_size = 8;
RMSE_info = "";
videos = ["cars" "flame"];
T = [3 5 7];
create_directory("../results");
for video = videos
    create_directory("../results/"+video);
    for t = T
        RMSE_info = process(video, t, dimensions, patch_size, seed, threshold, RMSE_info);
    end
end
file = fopen('../results/RMSE.txt', 'wt');
fprintf(file, RMSE_info);
fclose(file);
function RMSE_info = process(image_name, T, dimensions, patch_size, seed, threshold, RMSE_info)
    rng(seed);

    video_input = mmread(char("../results/"+ image_name +".avi"));
    for i=1:T
        X(:,:,i) = double(rgb2gray(video_input.frames(i).cdata));
    end
    % cropping
    % number_of_rows = dimensions(1);
    % number_of_columns = dimensions(2);
    % image_input = image_input(1:number_of_rows, 1:number_of_columns);
    [H,W,T] = size(X);
    code = randi([0 1], size(X));
    coded_snapshot = zeros(size(X));
    for i=1:T
        coded_snapshot(:,:,i) = X(:,:,i).*code(:,:,i);
    end
    coded_snapshot = sum(coded_snapshot, 3) + generate_gaussian_noise(dimensions, 0, 4);

    for i=1:T
        save_image(X(:,:,i), "../results/"+image_name+"/"+string(i)+".png")
    end
    save_image(coded_snapshot/T, "../results/"+image_name+"/"+"coded snapshot, T = "+string(T)+".png");
    save_image(coded_snapshot./(max(sum(code, 3), ones(dimensions))), "../results/"+image_name+"/"+"coded snapshot weighted averaging, T = "+string(T)+".png");

    % [image_reconstructed, RMSE] = patch_reconstruct(image_input, patch_size, measurement_matrix, sparse_basis_matrix, threshold);
    % create_directory("../results/"+basis);
    % save_image(image_reconstructed, "../results/" + basis + "/" + image_name + " " + comment + ".png");
    % RMSE_info = RMSE_info + "RMSE = " + string(RMSE) + " for " + image_name + " " + comment + " and " + basis + " basis " + newline;
    
    % RMSE = calculate_RMSE(image_input, image_with_gaussian_noise);
    % RMSE_info = RMSE_info + "RMSE = " + string(RMSE) + " for " + filename + newline;

    % RMSE_info = generate_result("reconstructed using compressive measurements, with noise", basis, image_name, image_with_gaussian_noise, patch_size, measurement_matrix, sparse_basis_matrix, threshold, RMSE_info);
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
            y = measurement_matrix * x(:);
            theta_hat = ista(y,  measurement_matrix*sparse_basis_matrix, threshold, 1);
            x_hat = sparse_basis_matrix * theta_hat;
            image_reconstructed(i:i+patch_size-1, j:j+patch_size-1) = image_reconstructed(i:i+patch_size-1, j:j+patch_size-1) + reshape(x_hat, [patch_size, patch_size]);
            number_of_overlapping_patches(i:i+patch_size-1, j:j+patch_size-1) = number_of_overlapping_patches(i:i+patch_size-1, j:j+patch_size-1) + ones(patch_size, patch_size);
        end
    end
    image_reconstructed = image_reconstructed ./ number_of_overlapping_patches;
    RMSE = calculate_RMSE(image_input, image_reconstructed)
end

function gaussian_noise_matrix = generate_gaussian_noise(size, mean, variance)
    gaussian_noise_matrix = sqrt(variance)*(mean + randn(size));
    % randn(size):samples 'size' elements from standar gaussian 
    % shifts these to N(mean,variance)
end

% Utility functions from my CS663 assignments
function RMSE = calculate_RMSE(image_original, image_reconstructed)
    RMSE = norm(image_reconstructed(:) - image_original(:)) / norm(image_original(:));
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