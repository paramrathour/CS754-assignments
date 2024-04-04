seed = 0;
%% Initialising Values
n=500; m=300;
L0_NORM = [5,10,15,20];
% L0_NORM = [5 10 ];
LAMBDA=[0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 30, 50, 100];

%% Do something

% Sensing Matrix one time generated
A = sensing_matrix_basis_matrix(m,n);
x_superset = sparsex(n,max(L0_NORM));
non_zero_indices=find(x_superset);
create_directory("../results")
for l0_norm=L0_NORM
    % select any l0_norm number of non-zero entries in the signal x_superset
    selected_indices = randsample(non_zero_indices,l0_norm);
    x_original =zeros(n,1); % select any l0_norm number of non-zero entries from the signal x_superset;
    x_original(selected_indices) = x_superset(selected_indices);
    % out of the m measurements,a random subset say 90% of the measurements are kept in Reconstructed Set R
    % and the remaining 10% are kept in Validation Set V
    indices_of_reconstruction_set = randperm(m,round(0.9*m));
    indices_of_validation_set = setdiff(1:m,indices_of_reconstruction_set);
    % I want to multiply x_original with with 1D DCT basis
    
    y_original = A*(x_original);
    noise_variance= 0.025*sum(abs(y_original))/m;
    noisy_y = y_original + sqrt(noise_variance)*randn(m,1);
    
    validation_error = zeros(size(LAMBDA));
    variance_error = zeros(size(LAMBDA));
    rmse_x = zeros(size(LAMBDA));
    rmse_x_e = zeros(size(LAMBDA));
    
    for lambda_index = 1:length(LAMBDA)
        lambda = LAMBDA(lambda_index);
        % the signal x is reconstructed using the measurements in R(and thus only the corresponding rows of sensing matrix)
        [reconstructed_x,status]=l1_ls(A(indices_of_reconstruction_set,:),noisy_y(indices_of_reconstruction_set),lambda);
        % the error in the reconstructed signal is calculated based on validation set
        % validation error corresponding to each lambda is calculated and stored
        % validation_error = norm(noisy_y(indices_of_validation_set)-A(indices_of_validation_set,:)*reconstructed_x)^2;
        validation_error(lambda_index)=(1.0/size(indices_of_validation_set,2))*norm(noisy_y(indices_of_validation_set)-A(indices_of_validation_set,:)*reconstructed_x)^2;
        disp(size(indices_of_validation_set,2));
        rmse_x(lambda_index) = norm(x_original-reconstructed_x)/norm(x_original);
        
        %solving part e
        % Pick lambda for which | ||y-sensing_matrix*reconstructed_x||^2-m*noise_variance | is minimum
        [reconstructed_x_e,status_e]=l1_ls(A,noisy_y,lambda);
        rmse_x_e(lambda_index) = norm(x_original-reconstructed_x_e)/norm(x_original);
        variance_error(lambda_index) = abs(norm(noisy_y-A*reconstructed_x_e)^2-m*noise_variance);
    end
    figure;
    semilogx(LAMBDA,validation_error,'o-');
    xlabel('Lambda(log scale)');
    ylabel('Validation Error');
    % title('Validation Error vs Log Lambda');
    title(['Validation Error vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    saveas(gcf, fullfile("../results", ['Validation_Error_Lambda_l0_norm_' num2str(l0_norm) '.png']));
    
    
    figure;
    semilogx(LAMBDA, rmse_x, '*-');
    xlabel('Lambda (log scale)');
    ylabel('RMSE');
    % title('RMSE vs. Logarithm of Lambda');
    title(['RMSE vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    saveas(gcf, fullfile("../results", ['RMSE_Lambda_l0_norm_' num2str(l0_norm) '.png']));
    
    % Plot rmse_x_e
    figure;
    semilogx(LAMBDA, rmse_x_e, 's-');
    xlabel('Lambda (log scale)');
    ylabel('RMSE_x_e');
    title(['RMSE_x_e vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    saveas(gcf, fullfile("../results", ['RMSE_x_Morozov_Lambda_l0_norm_' num2str(l0_norm) '.png']));
    
    % Plot variance_error
    figure;
    semilogx(LAMBDA, variance_error, '^-');
    xlabel('Lambda (log scale)');
    ylabel('Variance Error');
    % title(['Variance Error(Morozov's) vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    title(['Variance Error(Morozov) vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
   
    saveas(gcf, fullfile("../results", ['Variance_Error_Lambda_l0_norm_' num2str(l0_norm) '.png']));
    
    
end



%% Function to generate  sparse vector and sensing matrix
function x=sparsex(n,l0_norm)
x = zeros(n,1);
% draw non-zero entries of x at randomly chosen locations and let their values be drawn randomly from uniform (0,1000)
x(randperm(n,l0_norm)) = 1000*rand(l0_norm,1);
% randperm(n,k) returns a row vector containing k unique integers selected randomly from 1 to n
end
function M = sensing_matrix_basis_matrix(m,n)
% Sensing matrix should have entries -1 and +1 with equal probability
matrix =  2*(randi([0,1],m,n)-0.5);
matrix = matrix/sqrt(m);
M      = matrix*dctmtx(n);
end

%%
function [] = create_directory(path)
if ~exist(path)
    mkdir(path)
end
end