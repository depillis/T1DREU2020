%% TRANSLATING MCMC TUTORIAL ASTROSTATS
% Authors: R Code by C.A.L. Bailer-Jones, translated to MATLAB by Christina
%          Catlett, Maya Watanabe
% Date: 5/20/2020
% Description: R Code from Appendix C of 'Astrostats: Bayesian parameter estimation and
%              model comparison' translated to MATLAB


% Return a covariance natrix gievn a vector of the standard deviations and
% the global (i.e. scalr) correlation coefficient.
% This is inefficient, but matrix is very small and only done once.
% REFACTOR: make.covariance.matrix as covariance matrix
function covariancematrix = makecovmatrix(sampleSD, sampleCor)
    Len = length(sampleSD);
    covMatrix = zeros(Len, Len);
    if abs(sampleCor) > 1
        error('|sampleCor| > 1')
    end
    for i = 1:Len
        for j = 1:Len
            if i == j
                covMatrix(i,j) = sampleSD(i)^2;
            else
                covMatrix(i,j) = sampleCor*sampleSD(i)*sampleSD(j);
            end
        end
    end
    covariancematrix = covMatrix;
end
