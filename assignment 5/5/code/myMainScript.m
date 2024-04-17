seed = 0;
rng(seed);
n = 128;
number_of_samples = 10;
A = randn(n,n);
[U, ~] = qr(A);
c = 1;
alpha_values = 0:5;
m_values = [40 50 64 80 100 120];
RMSEs = zeros(length(alpha_values), length(m_values), length(number_of_samples));
for i = 1:length(alpha_values)
	alpha = alpha_values(i);
	eigenvalues = c*(1:n).^(-alpha);
	Sigma_x = U*diag(eigenvalues)*U';
	phi = generate_gaussian_noise(size(Sigma_x), 0, 1);
	samples = mvnrnd(zeros(n,1), Sigma_x, number_of_samples);
	for j = 1:length(m_values)
		m = m_values(j);
		variance = 1/m;
		phi_m = sqrt(variance)*phi(1:m, :);
		noise_vector = randn([m 1]);
		for k = 1:number_of_samples
			x = samples(k, :)';
			y = phi_m*x;
			sigma_noise = 0.01*mean(abs(y));
			noise = sigma_noise * noise_vector;
			y = y + noise;
			x_reconstructed = MAP(x, y, phi_m, Sigma_x, sigma_noise);
			RMSEs(i, j, k) = calculate_RMSE(x, x_reconstructed);
		end
	end
end
% RMSEs
RMSEs_average = mean(RMSEs, 3)

figure;
hold on;
for i = 1:length(alpha_values)
	plot(m_values, RMSEs_average(i, :), '-o');
end
hold off;
legend('alpha = ' + string(alpha_values));
saveas(gcf, "../../media/Q5 RMSEs.png");

figure_log = copyobj(gcf, 0);
set(gca(figure_log), 'YScale', 'log')
saveas(figure_log, "../../media/Q5 RMSEs log Y-axis.png");

function x_reconstructed = MAP(x, y, phi, Sigma, sigma_noise)
	m = size(phi, 1);
	x_reconstructed = inv(phi'*phi/sigma_noise^2 + inv(Sigma))*phi'*y/sigma_noise^2;
end

function RMSE = calculate_RMSE(signal_original, signal_reconstructed)
    RMSE = norm(signal_original - signal_reconstructed) / norm(signal_original);
end

function gaussian_noise_matrix = generate_gaussian_noise(size, mean, variance)
    gaussian_noise_matrix = sqrt(variance)*(mean + randn(size));
    % randn(size):samples 'size' elements from standar gaussian 
    % shifts these to N(mean,variance)
end