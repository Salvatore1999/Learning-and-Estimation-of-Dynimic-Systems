function [ gradJ ] = myLRCostFuncGrad(Phi,Y,Theta )
    F = mySigmoid(Phi*Theta);
    gradJ = Phi'*(F-Y);
end
