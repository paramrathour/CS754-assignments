threshold = 1e-4;
dimensions = [32 32];
% m = [900];
m = linspace(100,1000,10);
m_special = [500 700];
% k = [100];
k = [5 10 20 30 50 100 150 200];
k_special = [5 50 200];
seed = 9;

create_directory("../images")
process(m, m_special, k, k_special, dimensions, "omp", threshold, seed);
process(m, m_special, k, k_special, dimensions, "iht", threshold, seed);

function [] = process(m, m_special, k, k_special, dimensions, algorithm, threshold, seed)
    rng(seed);
    algorithm
    DCT_matrix = dctmtx(dimensions(1)*dimensions(2));
    column_indices = randperm(dimensions(1)*dimensions(2), k(end));
    coefficients = rand(k(end), 1);
    measurement_matrix = sign(randn(m(end), dimensions(1)*dimensions(2)));
    RMSE = zeros(length(k), length(m));
    create_directory("../images/" + algorithm)

    for i = 1:length(k)
        image_vector = DCT_matrix(:, column_indices(1:k(i))) * coefficients(1:k(i));
        filename = "../images/" + "k = " + string(k(i)) + ".png";
        maximum_brightness = max(image_vector, [], "all");
        save_image(image_vector, dimensions, filename, maximum_brightness);
        
        for j = 1:length(m)
            % disp([k(i) m(j)]) %
            phi = measurement_matrix(1:m(j), :);
            image_reconstructed = feval(algorithm, dimensions, phi*image_vector, k(i), phi*DCT_matrix, threshold);
            image_reconstructed = DCT_matrix * image_reconstructed;
            filename = "../images/" + algorithm + "/k = " + string(k(i)) + ", m = " + string(m(j)) + ".png";
            save_image(image_reconstructed, dimensions, filename, maximum_brightness);
            RMSE(i,j) = calculate_RMSE(image_vector, image_reconstructed, dimensions);
        end
    end
    fig1 = semilogy(m,RMSE(find(ismember(k, k_special)), :));
    legend("k = " + k_special)
    saveas(gcf, "../images/" + algorithm + "/plot " + "k" + ".png");
    fig2 = semilogy(k,RMSE(:, find(ismember(m, m_special))));
    legend("m = " + m_special)
    saveas(gcf, "../images/" + algorithm + "/plot " + "m" + ".png");
end

% Utility functions from my CS663 assignments
function RMSE = calculate_RMSE(image_original, image_reconstructed, dimensions)
    RMSE = norm(image_reconstructed - image_original) / norm(image_original);
end

function image_integer = image_double_to_integer(image_double)
    image_double(image_double < 0) = 0;
    image_double(image_double > 255)= 255;
    image_integer = uint8(image_double);
end

function [] = save_image(image, dimensions, filename, maximum_brightness)
    image = reshape(image, dimensions);
    image = image*255/maximum_brightness;
    image = image_double_to_integer(image);
    imwrite(image, filename);
end

function [] = create_directory(path)
    if ~exist(path)
        mkdir(path)
    end
end