function theta = ista(y, A, threshold, lamdba)
    theta = randn([size(A,2) 1]);
    % generates a random vector sampled from standard Gaussian with size
    % Number of columns of A by 1
    alpha = eigs(A'*A,1);
    % d=eigs(A,k) return k largest magnitude eigenvalues
    while(true)
        theta_next = wthresh(theta + 1/alpha*A'*(y-A*theta), 's', lamdba / (2*alpha));
        if norm(theta_next-theta) < threshold
            break
        end
        theta = theta_next;
    end
end