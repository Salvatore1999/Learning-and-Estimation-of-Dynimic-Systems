function [ MDL ] = myMDL( y , u, theta )

    M = length(y);
    p = length(theta);
    n = p/2;
    N = M - n;
    
    J = myCostFunc(y,u,theta);

    MDL = N*log(J)+2*p*log(N);

end

