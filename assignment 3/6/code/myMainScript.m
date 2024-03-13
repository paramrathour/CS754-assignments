seed = 0;
N = [50, 100, 500, 1000, 2000, 5000, 10000];
image_name = ["cryoem"];

rng(seed);
angles = 360*rand(max(N), 1);

for image = image_name
    for n = N
        image_input = imread("../results/"+ image +".png");
        image_input = double(image_input);
        projections = generateProjections(image_input, angles(1:n));
        image_output = reconstructImage(projections, angles(1:n));
        save_image(image_output, "../results/"+image+", N = "+string(n)+".png");
    end
end

function projections = generateProjections(image_input, angles)
    % returns a matrix where each column is a radon projections for corresponding angle
     projections = radon(image_input, angles);
end

function image_output = reconstructImage(projections, angles)
    % nearest neighbour algorithm
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