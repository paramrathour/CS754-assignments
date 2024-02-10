function matrix = haarWavelet(patch_size)
    [cA,cH,cV,cD] = dwt2(eye(patch_size*patch_size), 'haar');
    matrix = [cA cH; cV cD];
end