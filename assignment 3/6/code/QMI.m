% Note: This function is taken from my CS-663 assignments

function output = QMI(I1, I2)
    % we assume that the dimensions for both images are the same
    [height, width] = size(I1);

    p = JointHistogram(I1, I2) / (height * width);

    p_I1 = sum(p, 2); % sum of rows is for I1
    p_I2 = sum(p);    % sum of columns is for I2

    % now we have all the probability distributions

    output = 0;
    for i = 1:26
        for j = 1:26
            output = output + (p(i, j) - p_I1(i) * p_I2(j)) ^2;
        end
    end
end