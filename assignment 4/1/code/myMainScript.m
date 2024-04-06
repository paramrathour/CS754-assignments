seed = 0;
%% Initialising Values
n=500; m=300;
L0_NORM = [5,10,15,20];
% L0_NORM = [5,10];
LAMBDA=[0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 30, 50, 100];
% LAMBDA=[0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10];

%% Do something

% Sensing Matrix one time generated
A = sensing_matrix_basis_matrix(m,n);
x_superset = sparsex(n,max(L0_NORM));
non_zero_indices=find(x_superset);
% Open the text file for writing
fileID = fopen('results.txt', 'w');
validation_error_matrix = zeros(length(L0_NORM), length(LAMBDA));
rmse_x_matrix = zeros(length(L0_NORM), length(LAMBDA));
rmse_x_e_matrix = zeros(length(L0_NORM), length(LAMBDA));
variance_error_matrix = zeros(length(L0_NORM), length(LAMBDA));
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
    fprintf('lo_norm: %d\n', l0_norm);
    fprintf('Original y: %s\n', mat2str(y_original));
    fprintf('Noise variance: %.4f\n', noise_variance);
    fprintf('Noisy y: %s\n', mat2str(noisy_y));
    validation_error = zeros(size(LAMBDA));
    variance_error = zeros(size(LAMBDA));
    rmse_x = zeros(size(LAMBDA));
    rmse_x_e = zeros(size(LAMBDA));
    % following helper variables are defined to find out which lambda gave minimum error
    
    
    for lambda_index = 1:length(LAMBDA)
        lambda = LAMBDA(lambda_index);
        % the signal x is reconstructed using the measurements in R(and thus only the corresponding rows of sensing matrix)
        % [reconstructed_x]=l1_ls(A(indices_of_reconstruction_set,:),noisy_y(indices_of_reconstruction_set),lambda);
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
    % add validation_error, rmse_x, rmse_x_e, variance_error to the respective matrices
    validation_error_matrix(l0_norm==L0_NORM,:) = validation_error;
    rmse_x_matrix(l0_norm==L0_NORM,:) = rmse_x;
    rmse_x_e_matrix(l0_norm==L0_NORM,:) = rmse_x_e;
    variance_error_matrix(l0_norm==L0_NORM,:) = variance_error;

    [min_validation_error, arg_min_validation_error_index] = min(validation_error);
    [min_rmse_x, arg_min_rmse_x_index] = min(rmse_x);
    [min_rmse_x_e, arg_min_rmse_x_e_index] = min(rmse_x_e);
    [min_variance_error, arg_min_variance_error_index] = min(variance_error);
    
 
    
    % Check if the file was opened successfully
    if fileID == -1
        error('Unable to open file for writing');
    end
    % also write the l0_norm in text file
    fprintf(fileID, 'l0_norm: %d\n', l0_norm);
    % Write the results to the file
    fprintf(fileID, 'Minimum Validation Error: %.4f (Lambda: %.4f, Index: %d)\n', ...
        min_validation_error, LAMBDA(arg_min_validation_error_index), arg_min_validation_error_index);
    fprintf(fileID, 'Minimum RMSE_x: %.4f (Lambda: %.4f, Index: %d)\n', ...
        min_rmse_x, LAMBDA(arg_min_rmse_x_index), arg_min_rmse_x_index);
    fprintf(fileID, 'Minimum RMSE_x_e: %.4f (Lambda: %.4f, Index: %d)\n', ...
        min_rmse_x_e, LAMBDA(arg_min_rmse_x_e_index), arg_min_rmse_x_e_index);
    fprintf(fileID, 'Minimum Variance Error: %.4f (Lambda: %.4f, Index: %d)\n', ...
        min_variance_error, LAMBDA(arg_min_variance_error_index), arg_min_variance_error_index);
    
    % Close the file
    % fclose(fileID);
    fig1=figure;
    fig2=figure;
    fig3=figure;
    fig4=figure;

%{
     generateFigure(LAMBDA, validation_error, 'Validation Error', l0_norm, "../../media",fig1);
    generateFigure(LAMBDA, rmse_x, 'RMSE', l0_norm, "../../media",fig2);
    generateFigure(LAMBDA, rmse_x_e, 'RMSE-Morozov', l0_norm, "../../media",fig3);
    generateFigure(LAMBDA, variance_error, 'Morozov-variance-deviation', l0_norm, "../../media",fig4);
     
%}    
end
generateFigure(LAMBDA, validation_error_matrix, 'Validation Error', L0_NORM, "../../media");
generateFigure(LAMBDA, rmse_x_matrix, 'RMSE', L0_NORM, "../../media");
generateFigure(LAMBDA, rmse_x_e_matrix, 'RMSE-Morozov', L0_NORM, "../../media");
generateFigure(LAMBDA, variance_error_matrix, 'Morozov-variance-deviation', L0_NORM, "../../media");

fclose(fileID);
    
% modifying the generateFigure function such that it takes the matrices as input and plots the graph for all l0_norms
% function generateFigure(LAMBDA, y_axis_matrix,y_label, l0_norm, directory)
%     % Create directory if it doesn't exist
%     create_directory(directory)
%     % Create figure
%     figure(fig)
%     % fig = figure;
%     markers = {'o', '*', 's', '^', 'x', 'd', 'p', 'h'};
%     for i=1:size(y_axis_matrix,1)
%         marker_index = randi(length(markers));
%         % add legend corresponding to which l0_norm the curve belongs to
%         semilogx(LAMBDA, y_axis_matrix(i,:), [markers{marker_index} '-']); % Separate marker and line style
%         hold on


%     end
%     xlabel('Lambda (log scale)');
%     ylabel(y_label);
%     % title(['Validation Error vs. Logarithm of Lambda (l0_norm = ' num2str(l0_norm) ')']);


%{
 function generateFigure(LAMBDA,y_axis_matrix,y_label,L0_NORM,directory)
% Create directory if it doesn't exist L0_NORM
       create_directory(directory);
       figure;
       markers = {'o', '*', 's', '^', 'x', 'd', 'p', 'h'};
       for i=1:size(y_axis_matrix,1)
           marker_index= randi(length(markers));
           semilogx(LAMBDA,y_axis_matrix(i,:),[markers{marker_index} '-'],'DisplayName',['l0\_norm=' num2str(L0_NORM(i))]) ;
           hold on;
       end
       xlabel('Lambda(log scale)');
       ylabel(y_label);
       filename = fullfile(directory, ['Q1_' y_label '_l0_norm_' num2str(L0_NORM(i)) '.png']);
       saveas(gcf, filename);
end  
%}
function generateFigure(LAMBDA, y_axis_matrix, y_label, L0_NORM, directory)
    % Create directory if it doesn't exist L0_NORM
    create_directory(directory);
    figure;
    markers = {'o', '*', 's', '^', 'x', 'd', 'p', 'h'};
    for i = 1:size(y_axis_matrix, 1)
        marker_index = randi(length(markers));
        semilogx(LAMBDA, y_axis_matrix(i,:), [markers{marker_index} '-'], 'DisplayName', ['l0\_norm=' num2str(L0_NORM(i))]);
        hold on;
    end
    xlabel('Lambda (log scale)');
    ylabel(y_label);
    legend('Location', 'best'); % Add this line to display the legend
    filename = fullfile(directory, ['Q1_' y_label '_l0_norm_' num2str(L0_NORM) '.png']);
    saveas(gcf, filename);
end

%{
 function generateFigure(LAMBDA, y_axis_values,y_label, l0_norm, directory,fig)
% Create directory if it doesn't exist
% create_directory("../results")
create_directory(directory)
% Create figure
figure(fig)
% fig = figure;
markers = {'o', '*', 's', '^', 'x', 'd', 'p', 'h'};
marker_index = randi(length(markers));
semilogx(LAMBDA, y_axis_values, [markers{marker_index} '-']); % Separate marker and line style

xlabel('Lambda (log scale)');
ylabel(y_label);
% title(['Validation Error vs. Logarithm of Lambda (l0_norm = ' num2str(l0_norm) ')']);
legend(['$||\theta||_0$ = ' num2str(l0_norm)], 'Interpreter', 'latex');
% Save figure
filename = fullfile(directory, ['Q1_' y_label '_l0_norm_' num2str(l0_norm) '.png']);
saveas(fig, filename);
end 
%}



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

%{
    figure;
    semilogx(LAMBDA, rmse_x, '*-');
    xlabel('Lambda (log scale)');
    ylabel('RMSE');
    % title('RMSE vs. Logarithm of Lambda');
    title(['RMSE vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    saveas(gcf, fullfile("../results", ['RMSE_Lambda_l0_norm_' num2str(l0_norm) '.png']));
     
%}


%{
     % Plot rmse_x_e
    figure;
    semilogx(LAMBDA, rmse_x_e, 's-');
    xlabel('Lambda (log scale)');
    ylabel('RMSE_x_e');
    title(['RMSE_x_e vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    saveas(gcf, fullfile("../results", ['RMSE_x_Morozov_Lambda_l0_norm_' num2str(l0_norm) '.png']));
      lambda_rmse_x=0;
    lambda_validation=0;
    lambda_variance_morozorov=0;
    lambda_rmse_x_e=0;
    min_error_rmse_x=inf;
    min_error_validation=inf;
    min_error_variance_morozorov=inf;
    min_error_rmse_x_e=inf;
%}

%{
   
    % Plot variance_error
    figure;
    semilogx(LAMBDA, variance_error, '^-');
    xlabel('Lambda (log scale)');
    ylabel('Variance Error');
    % title(['Variance Error(Morozov's) vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    title(['Variance Error(Morozov) vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
     
%}

% saveas(gcf, fullfile("../results", ['Variance_Error_Lambda_l0_norm_' num2str(l0_norm) '.png']));

%{
     figure;
    semilogx(LAMBDA,validation_error,'o-');
    xlabel('Lambda(log scale)');
    ylabel('Validation Error');
    % title('Validation Error vs Log Lambda');
    title(['Validation Error vs. Logarithm of Lambda (l0\_norm = ' num2str(l0_norm) ')']);
    saveas(gcf, fullfile("../results", ['Validation_Error_Lambda_l0_norm_' num2str(l0_norm) '.png']));
     
%}