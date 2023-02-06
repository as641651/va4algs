function [res, time] = naive(X, M, y)
    tic;
    b = inv(transpose(X)*inv(M)*X)*transpose(X)*inv(M)*y;    
    time = toc;
    res = {b};
end