function [ J ] = myLRCostFunc(Phi,Y,Theta )
    N = length(Y);
    J = 0;
    for i=1:N
        phi = Phi(i,:);
        z = phi*Theta;
        J = J + Y(i)*log(mySigmoid(z))+(1-Y(i))*log(1-mySigmoid(z));
    end
end

