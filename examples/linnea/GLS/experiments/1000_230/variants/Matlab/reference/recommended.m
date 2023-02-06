function [res, time] = recommended(X, M, y)
    tic;
    b = ((transpose(X)*((M)\X))\transpose(X))*((M)\y);    
    time = toc;
    res = {b};
end