% Note: This function is taken from my CS-663 assignments

function array = JointHistogram(I1, I2)
    array = zeros(26, 26);
    
    for i = 1:numel(I1)
        intensity_x = I1(i);
        intensity_y = I2(i);

        bin_x = uint8(floor((intensity_x / 10) + 1));
        bin_y = uint8(floor((intensity_y / 10) + 1));

        array(bin_x, bin_y) = array(bin_x, bin_y) + 1;
    end
end