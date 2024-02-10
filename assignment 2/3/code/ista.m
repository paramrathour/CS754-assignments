function theta = ista(y, A, threshold)
    theta = randn([size(A,2) 1]);
    lamdba = 1;
    alpha = eigs(A'*A,1);
    while(true)
        theta_next = wthresh(theta + 1/alpha*A'*(y-A*theta), 's', lamdba / (2*alpha));
        if norm(theta_next-theta) < threshold
            break
        end
        theta = theta_next;
    end
end