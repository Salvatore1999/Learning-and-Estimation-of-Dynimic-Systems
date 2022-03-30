function [Test_Result] = myCrossCorrChi(eps_vector,u,N,m,alpha,order)
%MYCROSSCORRCHI Summary of this function goes here
%   It checks the cross-correlation assumption by using the Chi square
%   statistical test
    
    r_eps_vector = autocorr(eps_vector,'NumLags', m) * var(eps_vector);
    r_eps0 = r_eps_vector(1);
    % r_vector = r_vector(2:end);
    
    r_u_vector = autocorr(u(order+1:end),'NumLags', m) * var(u(order+1:end));
    r_u0 = r_u_vector(1);
    
    % From the theory, we know that the estimate of the covariance matrix
    % of the input is given by the transpose of the input Hankel Matrix
    % times the input Hankel matrix.
    % Noitce that now we don't consider just N-order samples as
    % always but N - (m+1-order)for the Hankel matrix.
    % It because, by looking at the theory we need to take the input until
    % the time sample (t-m) of our phi_u when we put it into the Hankel.
    Hu_m = myHank(u, m+1-order);
    Sigma_u_estimate =  1/N*(Hu_m'*Hu_m);
    % Notice that it is more or less the same to do, in matlab (uncomment):
%     Sigma_u_estimate = cov(Hu_m);
    
    
    % Once we have computed the estimate of the variance matrix of u, we
    % need to compute the cross-correlation between epsilon and u.
    % We can use the crosscorr matlab function but we must notice that this
    % crosscorr returns the normalized cross correlation, and so we need to
    % post multiply with the square root of the product of r_u(0) and
    % r_eps(0).
    r_epsu = crosscorr(eps_vector,u(order+1:N),'NumLags',m) * sqrt(r_eps0*r_u0);
    
    % For the same reason of the input Hankel matrix, we need to take into
    % account not all the samples of the cross correlation between epsilon
    % and u but only of the terms that are contained between the order of
    % the ARX model and the m (that is, in this case N/4).
    r_epsu = r_epsu(order:m);
    
    % Finally, compute the chi square value    
    x_chi = N*(r_epsu'*pinv(Sigma_u_estimate*N)*r_epsu)/r_eps0;
    
    alpha_opp = 1 - alpha;
    if x_chi <= chi2inv(alpha_opp,m)
        Test_Result = true;
    else
        Test_Result = false;
    end
    

end

