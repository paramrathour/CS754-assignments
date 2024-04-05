seed = 0;
epsilon = 0.1;
n = 200;
ranks = [1 2 3 5 10 15 20 25 50];
% ranks = [3 5 10 15];
fractions = [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
% fractions = [0.2 0.6];
lambdas = 10 .^ (2:4);
delta_k = 2;

rng(seed);
matrices = {rand([n n]), randn([n n])};
matrix_types = ["rand", "randn"];


for i = 1:length(matrices)
	X = matrices{i};
	matrix_type = matrix_types(i)
	RMSEs = zeros(length(ranks), length(fractions), length(lambdas));
	RMSEs = RMSEs(:);
	[I, J, Lambda] = ndgrid(1:length(ranks), 1:length(fractions), 1:length(lambdas));

	parfor k = 1:length(RMSEs) % replace "parfor" with "for" if you don't wish to install parallel computing toolbox
		r = ranks(I(k));
		f = fractions(J(k));
		lambda = lambdas(Lambda(k));
		RMSEs(k) = process_each_matrix(X, n, r, f, lambda, seed, epsilon, delta_k);
	end

	RMSEs = reshape(RMSEs, [length(ranks) length(fractions) length(lambdas)]);
	figure('Name', matrix_type);
	save("../../media/Q2 RMSEs "+ matrix_type + ".mat", "RMSEs");
	[optimal_RMSEs, optimal_lambda_indices] = min(RMSEs, [], 3);
	optimal_RMSEs
	optimal_lambdas = lambdas(optimal_lambda_indices)
	imagesc(fractions, ranks, optimal_RMSEs, [0 1]);
	set(gca, 'YDir', 'normal');
	set(gca, 'XTick', []);
	set(gca, 'YTick', []);
	xlabel('increasing $f$ values $\rightarrow$', 'Interpreter', 'latex');
	ylabel('increasing $r$ values $\rightarrow$', 'Interpreter', 'latex');
	colormap(flipud(gray));
	colorbar;
	saveas(gcf, "../../media/Q2 RMSEs "+ matrix_type + ".png");
end

function RMSE = process_each_matrix(X, n, r, f, lambda, seed, epsilon, delta_k)
	rng(seed);
	m = floor(f*n^2);
	mask = zeros([n n]);
	mask(randperm(n^2, m)) = 1;

	[U,S,V] = svd(X);
	X = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';
	M = X .* mask + generate_gaussian_noise(size(X), 0, (0.02*mean(abs(X),"all"))^2);
	X_reconstructed = SVT(M, mask, lambda, epsilon, delta_k);
	RMSE = calculate_RMSE(X, X_reconstructed);
end

function gaussian_noise_matrix = generate_gaussian_noise(size, mean, variance)
    gaussian_noise_matrix = sqrt(variance)*(mean + randn(size));
    % randn(size):samples 'size' elements from standar gaussian 
    % shifts these to N(mean,variance)
end

function RMSE = calculate_RMSE(image_original, image_reconstructed)
    RMSE = norm(image_original - image_reconstructed, "fro") / norm(image_original, "fro");
end