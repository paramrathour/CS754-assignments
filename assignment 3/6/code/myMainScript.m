seed = 0;
N = [50, 100, 500, 1000, 2000, 5000, 10000];
% N = [500]; % work with this for testing
RMSE_info = "";
image_name = ["cryoem"];
angles_rotation = 0:359;

rng(seed);
angles = 180*rand(max(N), 1); % would also have to compare with reverse if it was 360

for image = image_name
    for n = N
        image_input = imread("../results/"+ image +".png");
        image_input = double(image_input);
        projections = generateProjections(image_input, angles(1:n));
        image_output = reconstructImage(projections, angles(1:n), size(image_input)); % can use reconstructImageAliter for lower space complexity
        image_output = image_bound(image_output);
        [RMSE, image_reconstructed, angle_rotated] = calculate_RMSE_optimal(image_input, image_output, angles_rotation);
        if angle_rotated <= length(angles_rotation);
            reflected = "";
        else
            reflected = ", reflected";
        end
        angle = mod(angle_rotated, length(angles_rotation));
        if angle ~= 0
            angle = angle - 1;
            angle = 360 - angle;
        end
        save_image(image_output, "../results/"+image+", N = "+string(n)+".png");
        save_image(image_reconstructed, "../results/"+image+", N = "+string(n)+reflected+" and rotated by angle = "+string(angle)+".png");
        RMSE_info = RMSE_info + "RMSE = " + string(RMSE) + " for " + image_name + " and, N = "+string(n) + newline;
    end
end
file = fopen('../results/RMSE.txt', 'wt');
fprintf(file, RMSE_info);
fclose(file);

function projections = generateProjections(image_input, angles)
    % returns a matrix where each column is a radon projections for corresponding angle
     projections = radon(image_input, angles);
end

function image_output = reconstructImage(projections, angles, image_input_size)
    number_of_angles = length(angles);
    projections_sorted = zeros(size(projections));         % nearest neighbors columns one after another

    projections_sorted(:,1) = projections(:,1);
    projections(:,1)=[];
    
    % find nearest neighbor of ith ordered-projection -> that will become its i+1th ordered-projection 
    for i=1:number_of_angles-1
        norms = vecnorm(projections - projections_sorted(:,i)); % calculating columnwise norm
        [~,argminimum] = min(norms);                        % argminimum is nearest neighbor of ith ordered-projection
        projections_sorted(:,i+1) = projections(:,argminimum);
        projections(:,argminimum) = [];                         % delete the used column from minuend
    end
    uniform_angles = linspace(0, 180*(number_of_angles-1)/number_of_angles, number_of_angles);
    image_output = iradon(projections_sorted, uniform_angles, image_input_size(1));
    imshow(image_output, []);
end

function [RMSE, image_reconstructed, argmaximum] = calculate_RMSE_optimal(image_original, image_output, angles_rotation)
    % Note: This function is inspired from my CS-663 assignment 1
    image_output_flipped = fliplr(image_output);
    RMSE_values = [arrayfun(@(angle) calculate_RMSE(image_original, imrotate(image_output, angle, "bilinear", "crop")), angles_rotation) arrayfun(@(angle) calculate_RMSE(image_original, imrotate(image_output_flipped, angle, "bilinear", "crop")), angles_rotation)];
    [maximum, argmaximum] = min(RMSE_values);
    if argmaximum <= length(angles_rotation)
        image_reconstructed = imrotate(image_output, angles_rotation(argmaximum), "bilinear", "crop");
    else
        image_reconstructed = imrotate(image_output_flipped, angles_rotation(argmaximum-length(angles_rotation)), "bilinear", "crop");
    end
    RMSE = norm(image_reconstructed(:) - image_original(:)) / norm(image_original(:));
end

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

function image_output = reconstructImageAliter(projections, angles, image_input_size)
    number_of_angles = length(angles);
    uniform_angles = linspace(0, 180*(number_of_angles-1)/number_of_angles, number_of_angles);
    indices = 1:length(angles);
    indices_sorted = zeros(1, length(angles));
    
    indices_sorted(1) = 1;
    indices = indices(indices ~= 1);
    for i = 1:length(angles)-1
        s = projections(:, indices_sorted(i));
        Q = projections(:, indices);
        [~, argminimum]= min(vecnorm(Q - s));
        indices_sorted(i+1) = indices(argminimum);
        indices = indices(indices ~= indices(argminimum));
    end

    image_output = iradon(projections(:, indices_sorted), uniform_angles, image_input_size(1));
    imshow(image_output, []);
end