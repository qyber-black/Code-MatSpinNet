function [x,fitHist] = deoptim(objective,DIM,options)
% objective = objective function to be evaluated
% DIM  = number of control parameters
% options:
%   population     = size of population in each generation
%   lb,ub          = lower and upper bound on control parameters
%   mutation       = mutation parameter (see line 42)
%   crossover      = crossover parameter (see line 52)
%   maxGen         = max number of generations

% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
% SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
% SPDX-License-Identifier: AGPL-3.0-or-later

NPDIM = options.population;
bounds = [options.lb; options.ub];
F = options.mutation;
CR = options.crossover;
maxGen = options.maxGen;

k1 = diff(bounds);  % parameter range = upper bound - lower bound
k2 = bounds(1,:);     % lower bound
lower_bound = ones(NPDIM,1) * k2;
upper_bound = ones(NPDIM,1) * (k1 + k2);
diff_bound = ones(NPDIM,1) * k1;

fitness = zeros(1,NPDIM);
tempfit = zeros(1,NPDIM);

%Step One ------Initialize first NP and obj and error
NP = i4_to_hammersley_sequence(DIM, NPDIM, zeros(1,DIM), zeros(1,DIM), ones(1,DIM))' + rand(NPDIM,DIM)/nthroot(NPDIM,DIM);
NP = diff_bound .* NP + lower_bound;

for i=1:NPDIM
  fitness(i) = feval(objective,NP(i,:));
end

mutate = zeros(NPDIM,DIM);
Cross  = zeros(NPDIM,DIM);

%Step 2 -------While stopping criteria is not satisfied
%% is there no other termination condition, e.g., population has not
%% changed after N iterations, etc. ??
for G = 1:maxGen
  fitHist(G) = min(fitness);

  % Mutation Step -- could be vectorized but not crucial
  for i=1:NPDIM
    index = randi([1,NPDIM],3,1);
    % repeat until index vector has at least 3 distinct values
    % one of which is the index of the current individual i
    while ((length(unique(index))~=3)||(sum(index==i)~=0))
      index = randi([1,NPDIM],3,1);
    end
    %% this mutation is different from standard choices
    mutate(i,:) = NP(index(1),:) + F * (NP(index(2),:) - NP(index(3),:));
  end
  %check upper and lower bound here
  mutate = max(min(mutate,lower_bound), upper_bound);

  % CrossOver Step
  % Generate a trial vector for each target vector X(i,G)
  for i = 1:NPDIM
    Cross(i,:) = NP(i,:);
    jindex = randi([1,DIM],1,1);
    swap   = unique([jindex find(rand(1,DIM)<CR)]);
    Cross(i,swap) = mutate(i,swap);
  end

  % Selection Step
  for i=1:NPDIM
    tempfit(i) = feval(objective,Cross(i,:));
  end
  replace          = find(tempfit<fitness);
  fitness(replace) = tempfit(replace);
  for i = replace
    NP(i,:) = Cross(i,:);
  end

end

[fitHist(G),ind] = min(fitness);
x = NP(ind,:);
