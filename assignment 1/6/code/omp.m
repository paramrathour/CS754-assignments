function image_reconstructed = omp(dimensions, y, k, A, threshold)
    r = y;
    T = [];
    A = A ./ vecnorm(A);
    for i = 1:k;
        [~, j] = max(r'*A);
        T = union(T, j);
        theta = pinv(A(:, T))*y;
        r = y - A(:, T)*theta;
    end
    image_reconstructed = zeros(dimensions(1)*dimensions(2),1);
    image_reconstructed(T) = theta;
end