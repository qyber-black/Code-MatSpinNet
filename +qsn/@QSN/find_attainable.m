function [x,p,q,theta,OK,err,info] = find_attainable (obj,ispin,jspin,s)
  % [x,p,q,theta,OK,err,info] = find_attainable (r,ispin,jspin,s) - Find attainability parameters
  %
  % Uses GA optimisation to find parameters x such that
  %
  %  r.is_attainable (ispin,jspin,[x s])
  %
  % returns that the probability from ispin to jspin in a ring is attainable given a
  % scalar scaling parameter s.
  %

  % Call is_attainable to find dimension for x
  p = obj.is_attainable(ispin,jspin,1);
  Nbar = size(p,1);
  lb = ones(Nbar,1);
  range = 10000;

  if exist('gaoptimset')
    options = gaoptimset;
    options = gaoptimset(options,'Display', 'iter');
    options = gaoptimset(options,'PlotFcns', { @gaplotbestf });
    options = gaoptimset(options,'FitnessLimit', 0);
    options = gaoptimset(options,'TolFun', eps);
    options = gaoptimset(options,'Generations', 100);
    options = gaoptimset(options,'PopulationSize', 100);
    options = gaoptimset(options,'PopInitRange', [1;10000]);

    % Search
    [x,fval,exitflag,output,population,score] = ga(@(x) obj.opt_attainable(ispin,jspin,x,s),Nbar,[],[],[],[],lb,[],[],[],options);

  else

    % Defines the problem for the genetic algorithm
    problem.name = 'attainable';
    problem.minimization = true; % 1 for minimization and 0 for maximization
    % Tells the genetic algorithm which of your functions you want to use
    problem.generation_method = @attainable_generate_random;

    problem.mutation_method = @attainable_mutation_gaussian;
    problem.crossover_method = @attainable_crossover_polarized;
    problem.evaluate = @(x,y,z) obj.opt_attainable(ispin,jspin,x,s);
    problem.print = @attainable_print;

    % Defines the instance
    problem.n_var = Nbar; % Problem size
    problem.interval_center = zeros(1,problem.n_var) + range/2;
    problem.interval_size = ones(1,problem.n_var) * range/2;
    problem.reduction_rate = 1;

    population = GA(problem);

    x = population.best;
  end

  % Check results
  [p,q,theta,OK,err,info] = obj.is_attainable (ispin,jspin,[x s]);

  function [ x ] = attainable_generate_random( problem, ~ )
    x = rand(1,problem.n_var).*problem.interval_size-problem.interval_size/2+problem.interval_center;
  end

  function [ x ] = attainable_mutation_gaussian( x, problema, ~, ~, ~ )
    limite_inferior = problema.interval_center - 0.5*problema.interval_size;
    limite_superior = problema.interval_center + 0.5*problema.interval_size;
    x = x + 0.05*randn()*(limite_superior-limite_inferior);
  end

  function [ x ] = attainable_crossover_polarized(parents,problem,~,population)
    %CROSSOVER polarized with the parents
    x1 = population.ind{parents(1)};
    x2 = population.ind{parents(2)};

    limite_inferior = problem.interval_center - 0.5*problem.interval_size;
    limite_superior = problem.interval_center + 0.5*problem.interval_size;

    %gera indivíduo
    beta1 = rand();
    beta2 = rand();
    alpha = 1.4 * beta1 * beta2 - 0.2;
    x = alpha * x1 + (1-alpha) * x2;

    %aplica reflexão
    x = reflexao(x, limite_inferior, limite_superior);
  end

  function [ x ] = reflexao( x, limite_inferior, limite_superior )
    %REFLEXAO Aplica reflexão no ponto x de acordo com os limites
    if (sum(x<limite_inferior)>0)
      x = limite_inferior + abs(x-limite_inferior);
    end
    if (sum(x>limite_superior)>0)
        x = limite_superior - abs(x-limite_superior);
    end
  end

  function attainable_print( population, ~, settings, stats)

    %PRINT the best solution in the population
    figure(1);
    hold off;
    subplot(2,1,1);
    hold off;
    plot(population.best);
    title(['Current solution with fx = ',num2str(population.best_fx) ,' in generation ',num2str(stats.gen), ' (', num2str(population.t),' with no improvement)']);
    soma = population.ind{1};
    for i=2:settings.n_ind
        soma = soma + population.ind{i};
    end
    hold on;
    plot(soma./settings.n_ind,'r-');
    legend({'Best solution','Average Solution'});
    subplot(2,1,2);
    hist(population.fx(1:settings.n_ind));
    title('Variety of fx values');
    drawnow;
  end

end