function [Test_Result] = myCrossCorrAnderson(eps_vector,u,N,m,alpha,order)
    
    % It is simpler since to find the normalized cross-correlation, we
    % need only to use the function crosscorr.
    gamma_vector = crosscorr(eps_vector,u(order+1:end),'NumLags', m);
    gamma_vector = gamma_vector(2:end);
    
    bar_m = 0;
    alpha_opp = 1 - alpha/2;
    limit = norminv(alpha_opp)/sqrt(N);
    
    % Compute the Gaussian   
    
    for tau = 1:(m-1)
        if abs(gamma_vector(tau)) <= limit
        else
            bar_m = bar_m + 1;
        end
    end
    
    if (bar_m/(m-1)) <= alpha
        Test_Result = true;
    else
        Test_Result = false;
    end
    

end

