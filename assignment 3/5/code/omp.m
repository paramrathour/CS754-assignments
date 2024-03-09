function image_reconstructed = omp(y, A, threshold)
    m = size(A,1);
    r = y;
    T = [];
    i = 0;
    while (norm(r) > threshold & i < m)
        [~, j] = max(r'*A);
        T = union(T, j);
        theta = pinv(A(:, T))*y;
        r = y - A(:, T)*theta;
        i = i + 1;
    end
    image_reconstructed = zeros(size(A,2), 1);
    image_reconstructed(T) = theta;
end