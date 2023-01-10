function p = prob_binomial_log (s, X, D)
% LogProbability of H using bionomial distribition
%
%    p = prob_binomial_log (s, X, D)
%
%    s  - Spin measured
%    X  - Parameters [length-N coupling_strength_d]
%    D  - experimental data
%
%    p  - Log probability

% Ring with N nodes and uniform coupling d where we don't know N and/or d.
%
%   p_{ij}(t) = abs( sum_k <i|pi_k|j> e^{-i lambda_k t})
%
% where
%
%   <i|pi_k|j> = 2/N cos(2*pi k(i-j)/N)
%
% which simplifies to 2/N for all k if i=j. The eigenvalues should be
%
%  lambda_k(d) = 2d cos(2*pi k/N) for k=0,...,floor(N-1/2)+1.  So
%
%  p_{11}(t) = 2/N abs(sum_k e^{-i lambda_k(d) t})

% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019 Sophie M Shermer <lw1660@gmail.com>, Swansea University
% SPDX-License-Identifier: AGPL-3.0-or-later

% Probabilities from H for ring on 1,1 trace (ONLY - FIXME)

% FIXME: THIS IS NOT CORRECT
%lambda = 2 * X(2) * cos (2 * pi * [0:1:floor(X(1)-.5)+1] / X(1));
%theta = arrayfun (@(t) 2 * abs(sum(exp(-i * lambda * t))) / X(1), D(1,:));

N = qsn.QSN('ring', ones(1,X(1)) * X(2));
theta = N.trace (s, s, D(1,:));

% Binomial probability of seeing data
w = arrayfun(@(n,k)nchoosek(n,k),D(2,:)+D(3,:),D(2,:));
lw = log(w);
lw(w <= eps(0.0)) = 0;

lt1 = log(theta);
lt1(theta <= eps(0.0)) = 0;
lt2 = log(1-theta);
lt2(1-theta <= eps(0.0)) = 0;

p = sum(lw + lt1.*D(2,:) + lt2.*D(3,:));

p = max(p -log(eps(0.0)), 0);
