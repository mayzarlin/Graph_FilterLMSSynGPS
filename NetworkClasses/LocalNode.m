classdef LocalNode < handle
    properties
        L  %% Local Laplacian
        S   %% Permuted Local EigenVector Matrix
        D   %% Permuted Local EigenValue Matrix
        EstS 
        C   %% Permuted Covariance Matrix 
        E   %% Permutation Matrix
    end
    methods 
        function Local = LocalNode(n)
            Local.L     = zeros(n);
            Local.S     = zeros(n);
            Local.D     = zeros(n);
            Local.EstS  = zeros(n);
            Local.C     = zeros(n);
            Local.E     = zeros(n);
        end
    end
end