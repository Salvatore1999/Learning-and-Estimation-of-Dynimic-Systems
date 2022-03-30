
% Hankel matrix generator
% function H = hank(X,n)
% s --> data vector
% n --> order

function H = myHank(s,n)
    M=length(s);
    N = M - n;
    H=zeros(N,n); %Init the matrix H

    for i = 1:(N)
        for j = 1:n
            H(i,j) = s(n-j+i);
        end
    end
end


