function matrix = dct2D(n)
    matrix = kron(dctmtx(n)',dctmtx(n)');
end