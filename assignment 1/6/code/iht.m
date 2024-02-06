function x = iht(dimensions, y, k, A, threshold)
    [U, S, V] = svd(A);
    scale_factor = 1/(2*max(diag(S)));
    phi = A * scale_factor;
    y = y * scale_factor;
    x = zeros(dimensions(1)*dimensions(2), 1);
    while(true)
        x_next = H(x+phi'*(y-phi*x), k);
        if (norm(x_next-x) < threshold)
            break;
        end
        x = x_next;
    end
end

function x = H(x, s)
    [~, indices] = sort(x, 'descend');
    x(~ismember(1:length(x), indices(1:s))) = 0;
end