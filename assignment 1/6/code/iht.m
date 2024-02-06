function x = iht(dimensions, y, k, phi, threshold)
    x = zeros(dimensions(1)*dimensions(2), 1);
    while(norm(y-phi*x) > threshold)
        x = H(x+phi'*(y-phi*x), k);
    end
end

function x = H(x, s)
    [~, indices] = sort(x, 'descend');
    ismember(1:length(x), indices(1:s));
    x(~ismember(1:length(x), indices(1:s))) = 0;
end