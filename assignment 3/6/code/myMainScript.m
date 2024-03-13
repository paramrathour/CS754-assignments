seed = 0;
N = [50, 100, 500, 1000, 2000, 5000, 10000];
N = [50]; % work with this for testing
RMSE_info = "";
image_name = ["cryoem"];
angles_rotation = 0:359;

rng(seed);
angles = 360*rand(max(N), 1);

for image = image_name
    for n = N
        image_input = imread("../results/"+ image +".png");
        image_input = double(image_input);
        projections = generateProjections(image_input, angles(1:n));
        image_output = reconstructImage(projections, angles(1:n), image_input);
        [RMSE, image_reconstructed, angle_rotated] = calculate_RMSE(image_input, image_output, angles_rotation);
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

function image_output = reconstructImage(projections, angles, image_input)
    % nearest neighbour algorithm
end

function [RMSE, image_reconstructed, argmaximum] = calculate_RMSE(image_original, image_output, angles_rotation)
    % Note: This function is inspired from my CS-663 assignment 1
    image_output_flipped = fliplr(image_output);
    QMI_values = [arrayfun(@(angle) QMI(image_original, imrotate(image_output, angle, "bilinear", "crop")), angles_rotation) arrayfun(@(angle) QMI(image_original, imrotate(image_output_flipped, angle, "bilinear", "crop")), angles_rotation)];
    [maximum, argmaximum] = max(QMI_values);
    if argmaximum <= length(angles_rotation)
        image_reconstructed = imrotate(image_output, angles_rotation(argmaximum), "bilinear", "crop");
    else
        image_reconstructed = imrotate(image_output_flipped, angles_rotation(argmaximum-length(angles_rotation)), "bilinear", "crop");
    end
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