function [Test_Result] = myWhitenessAnderson(eps_vector,N,m,alpha)
     
    % Now we are interested to the normalized autocorrelation vector.
    % in the slides we called it gamma, for notation I will use
    % gamma_vector
    
    gamma_vector = autocorr(eps_vector,'NumLags', m);
    gamma_vector = gamma_vector(2:end);
    bar_m = 0;
    alpha_opp = 1 - alpha/2;
    limit = norminv(alpha_opp)/sqrt(N);
    
    for tau = 1:(m-1)
        if abs(gamma_vector(tau)) <= limit
        else
            bar_m = bar_m + 1;
        end
    end
    
    if (bar_m /(m-1)) <= alpha
        Test_Result = true;
    else
        Test_Result = false;
    end
    
end

