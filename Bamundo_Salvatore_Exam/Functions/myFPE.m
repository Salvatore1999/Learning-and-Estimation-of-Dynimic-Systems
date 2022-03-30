function [ FPE ] = myFPE( y, u, theta )

    M = length(y);
    p = length(theta);
    n = p/2;
    N = M - n;
    
    J=myCostFunc(y,u,theta);

    FPE=(N+p)*J/(N-p);

end

