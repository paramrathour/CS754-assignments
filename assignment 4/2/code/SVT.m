function X_reconstructed = SVT(M, mask, lambda, epsilon, delta_k, maximum_iterations)
	Y = zeros(size(M));
	i = 0;
	while (i < maximum_iterations)
		X_reconstructed = svd_threshold(Y, lambda);
		Y = Y + delta_k * mask .* (M - X_reconstructed);
		if (norm(mask .* (X_reconstructed - M), "fro")/norm(mask .* M, "fro") <= epsilon) 
			break;
		end
		i = i + 1;
	end
end

function Y_soft = svd_threshold(Y, lambda)
	[U,S,V] = svd(Y);
	n = length(diag(S));
	S(1:n, 1:n) = diag(max(0, diag(S) - lambda));
	Y_soft = U * S * V';
end