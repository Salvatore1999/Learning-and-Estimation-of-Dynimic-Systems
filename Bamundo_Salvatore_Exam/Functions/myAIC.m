function [ AIC ] = myAIC( y , u, theta )

    M = length(y);
    p = length(theta);
    n = p/2;
    N = M - n;
    
    J = myCostFunc(y,u,theta);

    AIC = N*log(J)+2*p;

end

