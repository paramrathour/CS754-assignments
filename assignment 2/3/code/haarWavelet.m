function matrix = haarWavelet(n)
    I = eye(n*n);
    matrix = zeros([n*n n*n]);
    for i = 1:n*n
        vector = I(:, i);
        temp = dwt2Haar(reshape(vector, [n n]));
        matrix(:, i) = temp(:);
    end
end

function matrix = dwt2Haar(image)
    [cA,cH,cV,cD] = dwt2(image, 'haar');
    matrix = [cA cH; cV cD];
end