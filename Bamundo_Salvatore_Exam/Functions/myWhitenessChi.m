function [Test_Result]=myWhitenessChi(eps_vector,N,m,alpha)

    % remember that autocorr(-) gives us the normalized autocorrelation, but we
    % need of the normal autocorrelation for this step.
    r_vector = autocorr(eps_vector,'NumLags', m) * var(eps_vector);
    
    % From the theory we know that the vector of r goes from r(1) that in 
    % matlab is r(2)up to r(m) and so r(0) that in matlab is r(1) is not
    % contained. 
    r0 = r_vector(1);
    r_vector = r_vector(2:end);
    
    % slide 4/9 of "model_vadilation" formula:
    x_chi = N* (r_vector'*r_vector)/(r0^2);
    
    % To make the statistical test, I used "chi2inv(p,nu) that returns
    % the inverse cumulative distribution function (icdf) of the chi-square
    % distribution with degrees of freedom nu, evaluated at the probability
    % values in p".
    % https://it.mathworks.com/help/stats/chi2inv.html
    
    % Of course, our "p" is 1-alpha since we need to find the the internal
    % values of the chi2 distribution.
    
    % Mathworks example: chi2inv(95,10) = 18.3, then it means that
    % "If you generate random numbers from this chi-square distribution,
    % you would observe numbers greater than 18.3 only 5% of the time".  
    
    alpha_opp = 1 - alpha;
    if x_chi <= chi2inv(alpha_opp,m)
        Test_Result = true;
    else
        Test_Result = false;
    end

end