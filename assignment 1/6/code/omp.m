function image_reconstructed = omp(dimensions, y, k, measurement_matrix, DCT_matrix)
    r = y;
    T = [];
    A = measurement_matrix;
    A = A ./ vecnorm(A);
    for i = 1:k;
        [~, j] = max(r'*A);
        T = union(T, j);
        theta = pinv(A(:, T))*y;
        r = y - A(:, T)*theta;
    end
    image_reconstructed = zeros(dimensions(1)*dimensions(2),1);
    image_reconstructed(T) = theta;
    image_reconstructed = DCT_matrix * image_reconstructed;
end