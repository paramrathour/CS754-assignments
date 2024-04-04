seed = 0;
epsilon = 0.1;
n = 200;
ranks = [1 2 3 5 10 15 20 25 50];
% ranks = [3 5 10 15 20];
fractions = [0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8];
% fractions = [0.2 0.6];
lambdas = 10 .^ (0:3);
% lambdas = [100 1000];
delta_k = 2;

RMSEs = zeros(length(ranks), length(fractions), length(lambdas));
RMSEs = RMSEs(:);
[I, J, Lambda] = ndgrid(1:length(ranks), 1:length(fractions), 1:length(lambdas));

rng(seed);
X = rand([n n]);

parfor k = 1:length(RMSEs)
	r = ranks(I(k));
	f = fractions(J(k));
	lambda = lambdas(Lambda(k));
	RMSEs(k) = process(X, n, r, f, lambda, seed, epsilon, delta_k);
end

RMSEs = reshape(RMSEs, [length(ranks) length(fractions) length(lambdas)]);
save("RMSEs.mat", "RMSEs");
[optimal_RMSEs, optimal_lambda_indices] = min(RMSEs, [], 3);
optimal_RMSEs
optimal_lambdas = lambdas(optimal_lambda_indices)
imagesc(fractions, ranks, optimal_RMSEs);
% imagesc('XData', fractions, 'YData', ranks, 'CData', optimal_RMSEs);
set(gca, 'YDir', 'normal');
set(gca, 'XTick', []);
set(gca, 'YTick', []);
xlabel('$f$ values', 'Interpreter', 'latex');
ylabel('$r$ values', 'Interpreter', 'latex');
colormap(flipud(gray));
colorbar;
saveas(gcf, 'RMSEs.pdf');

function RMSE = process(X, n, r, f, lambda, seed, epsilon, delta_k)
	rng(seed);
	m = f*n^2;
	mask = zeros(n);
	mask(randperm(n^2, m)) = 1;

	[U,S,V] = svd(X);
	X = U(:, 1:r) * S(1:r, 1:r) * V(:, 1:r)';
	M = X .* mask + generate_gaussian_noise(size(X), 0, (0.02*mean(abs(X),"all"))^2);
	X_reconstructed = SVT(M, lambda, mask, epsilon, delta_k);
	RMSE = calculate_RMSE(X, X_reconstructed);
end

function gaussian_noise_matrix = generate_gaussian_noise(size, mean, variance)
    gaussian_noise_matrix = sqrt(variance)*(mean + randn(size));
    % randn(size):samples 'size' elements from standar gaussian 
    % shifts these to N(mean,variance)
end

function RMSE = calculate_RMSE(image_original, image_reconstructed)
    RMSE = norm(image_reconstructed(:) - image_original(:)) / norm(image_original(:));
end