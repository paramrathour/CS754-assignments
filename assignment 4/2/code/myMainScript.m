seed = 0;
epsilon = 1e-2;
n = 200;
ranks = [1 2 3 5 10 15 20 25 50];
fractions = [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
lambdas = [10 50 100 200];
delta_k = 2;

fraction_validation = 0.9;

rng(seed);
matrices = {randn([n n]), rand([n n])};
matrix_types = ["randn", "rand"];


for i = 1:length(matrices)
	X = matrices{i};
	[U,S,V] = svd(X);
	matrix_type = matrix_types(i)
	RMSEs = zeros(length(ranks), length(fractions), length(lambdas));
	validation_errors = zeros(length(ranks), length(fractions), length(lambdas));
	[I, J, Lambda] = ndgrid(1:length(ranks), 1:length(fractions), 1:length(lambdas));
	RMSEs = RMSEs(:);
	validation_errors = validation_errors(:);

	% the code is designed to be parallelisable
	% replace "for" with "parfor" if you want faster computation (needs parallel computing toolbox)
	parfor k = 1:length(RMSEs)
		r = ranks(I(k));
		f = fractions(J(k));
		lambda = lambdas(Lambda(k));
		X_r = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';
		[validation_errors(k), RMSEs(k)] = process_each_matrix(seed, X_r, f, fraction_validation, lambda, epsilon, delta_k);
	end

	RMSEs = reshape(RMSEs, [length(ranks) length(fractions) length(lambdas)]);
	validation_errors = reshape(validation_errors, [length(ranks) length(fractions) length(lambdas)]);
	figure('Name', matrix_type);
	save("../../media/Q2 "+ matrix_type + " RMSEs.mat", "RMSEs");
	save("../../media/Q2 "+ matrix_type + " Validation Errors.mat", "validation_errors");
	[optimal_validation_errors, optimal_lambda_indices] = min(validation_errors, [], 3);
	optimal_validation_errors
	optimal_RMSEs = RMSEs(sub2ind(size(RMSEs), repmat((1:length(ranks))', [1, length(fractions)]), repmat((1:length(fractions)), [length(ranks), 1]), optimal_lambda_indices))
	optimal_lambdas = lambdas(optimal_lambda_indices)
	imagesc(fractions, ranks, optimal_RMSEs, [0 1]);
	set(gca, 'YDir', 'normal');
	set(gca, 'XTick', []);
	set(gca, 'YTick', []);
	xlabel('increasing $f$ values $\longrightarrow$', 'Interpreter', 'latex');
	ylabel('increasing $r$ values $\longrightarrow$', 'Interpreter', 'latex');
	colormap(flipud(gray));
	colorbar;
	saveas(gcf, "../../media/Q2 "+ matrix_type + " RMSEs.png");
end

function [validation_error, RMSE] = process_each_matrix(seed, X, f, fraction_validation, lambda, epsilon, delta_k)
	rng(seed);
	[n_1, n_2] = size(X);
	m = floor(f*n_1*n_2);
	v = floor(fraction_validation*m);
	mask_measurements = zeros(size(X));
	mask_validation = zeros(size(X));
	mask_reconstruction = zeros(size(X));
	indices_measurement = randperm(n_1*n_2, m);
	indices_validation = indices_measurement(1:v);
	indices_reconstruction = indices_measurement(v+1:end);
	mask_measurements(indices_measurement) = 1;
	mask_validation(indices_validation) = 1;
	mask_reconstruction(indices_reconstruction) = 1;

	N = generate_gaussian_noise(size(X), 0, (0.02*mean(abs(X), "all"))^2);
	M_measurements = X .* mask_measurements + N;
	M_validation = X .* mask_validation + N;
	X_reconstructed_with_validation = SVT(M_validation, mask_validation, lambda, epsilon, delta_k);
	X_reconstructed_with_measurements = SVT(M_measurements, mask_measurements, lambda, epsilon, delta_k);

	validation_error = mean(mask_reconstruction .* (X - X_reconstructed_with_validation).^2, "all");
	RMSE = calculate_RMSE(X, X_reconstructed_with_measurements);
end

% Utility functions from my previous CS754 assignments
function gaussian_noise_matrix = generate_gaussian_noise(size, mean, variance)
    gaussian_noise_matrix = sqrt(variance) * (mean + randn(size));
    % randn(size):samples 'size' elements from standar gaussian 
    % shifts these to N(mean,variance)
end

function RMSE = calculate_RMSE(image_original, image_reconstructed)
    RMSE = norm(image_original - image_reconstructed, "fro") / norm(image_original, "fro");
end