function [ J ] = myCostFunc(y,u,theta)
    
    M = length(y);
    p = length(theta);
    n = p/2;
    N = M - n;    

    Hu = myHank(u,n); Hy = myHank(y,n); 

    H = [-Hy, Hu];
    
    J=1/N*(y(n+1:end)-H*theta)'*(y(n+1:end)-H*theta);

end

