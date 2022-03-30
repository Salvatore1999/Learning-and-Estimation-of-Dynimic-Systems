function [ HessJ ] = myLRCostFuncHessian(Phi,Y,Theta )
    F = mySigmoid(Phi*Theta);
    W = diag(F.*(1-F));
    HessJ = Phi'*W*Phi;
end

