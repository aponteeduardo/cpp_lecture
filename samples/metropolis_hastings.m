function [post] = metropolis_hastings(y, u, ptheta, htheta)
%% Skeleton for implementing the metropolis hasting algorithm
%
% Input
%   y           Experimental data
%   u           Experimental design
%   ptheta      Priors of the parameters
%   htheta      Parameters of the Metropolis hastings algorithm
%
% Output
%
% aponteeduardo@gmail.com
% copyright (C) 2015
%

% This gives different random number everytime
rng('shuffle')

% Number of burn in iteration
nburnin = htheta.nburnin;

% Number of iterations use to sample
niter = htheta.niter;

% Compute the square root of the covariance for the proposal distribution
htheta.csigma = chol(htheta.sigma);

% Samples from the posterior distribution
post = struct();
post.theta = zeros(numel(otheta), niter);

otheta = htheta.csigma * randn(size(htheta, 1), 1);
ollh = llh(y, u, otheta, ptheta);
olpp = lpp(u, theta, ptheta);

for i = 1: nburnin + niter
    % Propose a new sample
    ntheta = otheta + htheta.csigma * randn(size(otheta));
    
    % Compute log likelihood conditioned on a proposed sample
    nllh = llh(y, u, ntheta, ptheta);
    % Compute log prior probability
    nlpp = lpp(u, ntheta, ptheta);

    % Rejection criteriai in log space
    preject = max(0, nllh + nlpp - ollh - olpp);

    % If the samples is accepted keep it 
    if preject > log(rand)
        otheta = ntheta;
        ollh = nllh;
        olpp = nlpp;
    end

    % Store your samples
    if i > nburnin
        post.theta(:, i - nburnin) = ntheta;
    end

end


end % Metropolis hastings

function [l] = llh(y, u, theta, ptheta)
% A simple model of learning a policy

alpha = theta(1);
alpha = arctan(alpha)/pi + 0.5 % Squeeze to the [0, 1] interval

beta = theta(2);
beta = exp(beta); % Make it always positive.

vA = ptheta.vA0;
vB = ptheta.vB0;

l = 0;

for i = 1:numel(y)
    if y(i) == 0 % Case where decision is A
        % Likelihood of a decision
        l = l + beta * vA - log(exp(beta * vA) + exp(beta * vB));
        vA = vA + alpha(u(i) - vA);
    elseif y(i) == 1 % Case where decision is B
        % Likelihood of a decision
        l = l + beta * vB - log(exp(beta * vA) + exp(beta * vB));
        vB = vB + alpha(u(i) - vB);
    end
end

end % llh

function [p] = lpp(u, theta, ptheta)
% Prior probability of the parameters

% We asssume that alpha is beta distributed
alpha = theta(1);
alpha = arctan(alpha)/pi + 0.5 % Squeeze to the [0, 1] interval

p = log(betapdf(alpha, ptheta.alpha_A, ptheta.alpha_B));

% Beta is log normal distributed
p = p + log(normpdf(beta, ptheta.beta_mu, ptheta.beta_sigma));

end % lpp
