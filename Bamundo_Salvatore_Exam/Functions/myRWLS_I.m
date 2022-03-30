function [theta]=myRWLS_I(y,u,n,lambda)

    % In this case, length(y) is what I called before M
    M = length(y);
    N = M - n;
    % theta is the number of parameters that in ARX case is 2*n where n is
    % the order, the initialization is simply zeros
    theta0=zeros(2*n,1);
    phi1 = zeros(2*n,1);
    phit = zeros(2*n,1);
    
    % I create he firt regressor phi(1)=[-y(0)...-y(1-n) u(0)...u(1-n)]^T
    for i=1:n
        phi1(i) = -y(n-i+1);
        phi1(i+n) = u(n-i+1);
    end
    
    % Update S(1)
    S1=(phi1*phi1');
    % Compute K(1)
    K1 =  pinv(S1)*phi1;
    % Compute the error epsilon(1)
    epsilon1 = y(n+1)-phi1'*theta0;
    % update theta1
    theta1 = theta0 + K1*epsilon1;
    
    S_before=S1;
    theta_before = theta1;
    
    for i=2:N
        % phit is phi(t)
        for j=1:n
            phit(j) = -y(n-j+i);
            phit(j+n) = u(n-j+i);
        end
        
        % Update S(t)
        St = lambda*S_before+phit*phit';
        % Compute K(t)
        Kt =  pinv(St)*phit;
        % Compute the error epsilon(t)
        epsilont = y(n+i)-phit'*theta_before;
        % update theta1
        thetat= theta_before + Kt*epsilont;
        
        theta_before = thetat;
        S_before=St;
        
    end
    theta=thetat;
end