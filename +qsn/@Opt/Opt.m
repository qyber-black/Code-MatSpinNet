classdef (Sealed = true) Opt < handle
  % Opt   Optimiser class
  %
  % This class represents an optimiser object.
  %
  % OptSwitch Properties:
  %   alg           - Optimisation procedure
  %   problem       - Optimiser options
  %
  % OptSwitch Methods:
  %   obj = Opt (e, opt)                    - Optimiser c'tor
  %   [x,err,exec_time,exit_flag,output] = obj . solve (eval_func, x0)
  %                                         - Solve problem
  %
  properties (SetAccess = private)
    alg       = 'undef';                % Optimisation procedure
    tol       = double(eps('single'));  % Tolerance
    display   = 'off';                  % Display
  end
  methods
    function obj = Opt (a, e, dbg)
      % obj = Opt (a, e, opt) - Optimiser c'tor
      %
      %   a         - Optimisation procedure:
      %                'search'   - fminsearch
      %                'fmin'     - fminunc/fmincon
      %                'fmingrad' - fminunc/fmincon with gradient
      %                'ga'       - genetic algorithm
      %                'de'       - differential evolution
      %                'pattern'  - patternsearch
      %   e         - Precision (optional)
      %   dbg       - vlaue for Display option (default 'off')
      %   obj       - Optimiser object
      %
      switch a
        case 'search'
        case 'fmin'
        case 'fmingrad'
        case 'ga'
        case 'de'
        case 'pattern'
        otherwise
          error ('Illegal optimisation procedure selected');
      end
      obj.alg = a;
      if exist('e', 'var')
        obj.tol = e;
      end
      if exist('dbg', 'var')
        obj.display = dbg;
      end
    end
    function disp (obj)
      % obj.disp () - Display Opt object
      %
      %   obj - Opt object
      %
      fprintf (1, 'Optimiser: %s\nTolerance: %g\nDisplay: %s\n', obj.alg, obj.tol, obj.display);
    end
    function [x,err,exec_time,exit_flag,output] = solve (obj, eval_func, x0, cons)
      % [x,err,exec_time,exit_flag,output] = obj . solve ()
      %
      % Try to find solution by switching controls on and off.
      %
      %   eval_func - Target function
      %   x0        - Initial value
      %   obj       - Opt object
      %   x         - Switching times
      %   err       - Error
      %   exec_time - Execution time
      %   exit_flag - Optimisation exit flag status indicator
      %   output    - Info about optimisation
      %   x0        - Initial switching times
      %   cons      - constriants
      %
      % FIXME: better way to report results?
      if exist('cons','var')
        problem.Aineq   = cons.Aineq;
        problem.bineq   = cons.bineq;
        problem.Aeq     = cons.Aeq;
        problem.beq     = cons.beq;
        problem.lb      = cons.lb;
        problem.ub      = cons.ub;
        problem.nonlcon = cons.nonlcon;
        switch obj.alg
          case 'search'
            error('Search optimiser does not accept constraints');
          case 'fmin'
            problem.solver      = 'fmincon';
            opts                = optimoptions('fmincon');
            opts.Algorithm      = 'interior-point';
            opts.GradObj        = 'off';
            opts.ConstraintTolerance = obj.tol;
          case 'fmingrad'
            problem.solver      = 'fmincon';
            opts                = optimoptions('fmincon');
            opts.Algorithm      = 'trust-region-reflective';
            opts.GradObj        = 'on';
            opts.ConstraintTolerance = obj.tol;
          case 'ga'
            problem.solver      = 'ga';
            opts.PopInitRange   = [0;10];
            opts.PopulationSize = 200;
            opts.Generations    = 1000;
            opts.TolCon     = obj.tol;
          case 'de'
            problem.solver    = 'de';
            opts.population   = 100;
            opts.mutation     = 0.5;
            opts.crossover    = 0.5;
            opts.maxGen       = 100;
            opts.TolCon     = obj.tol;
          case 'pattern'
            problem.solver      = 'patternsearch';
            opts.Cache          = 'on';
            opts.CacheTol       = obj.tol;
            opts.TolMesh        = obj.tol;
            opts.TolCon     = obj.tol;
        end
      else
        problem.Aineq   = [ -eye(length(x0)) ];
        problem.bineq   = [ zeros(size(x0)) ];
        problem.Aeq     = [];
        problem.beq     = [];
        problem.lb      = [];
        problem.ub      = [];
        problem.nonlcon = [];
        switch obj.alg
          case 'search'
            problem.solver      = 'fminsearch';
          case 'fmin'
            problem.solver      = 'fminunc';
            opts                = optimoptions('fminunc');
%            opts                = optimoptions('fmincon');
            opts.Algorithm      = 'quasi-newton';
%            opts.LargeScale     = 'off';
            opts.GradObj        = 'off';
          case 'fmingrad'
            problem.solver      = 'fminunc';
            opts                = optimoptions('fmincon');
            opts.Algorithm      = 'active-set';
            opts.GradObj        = 'on';
          case 'ga'
            problem.solver      = 'ga';
            opts.PopInitRange   = [0;10];
            opts.PopulationSize = 200;
            opts.Generations    = 1000;
          case 'de'
            problem.solver    = 'de';
            opts.population   = 100;
            opts.mutation     = 0.5;
            opts.crossover    = 0.5;
            opts.maxGen       = 100;
          case 'pattern'
            problem.solver      = 'patternsearch';
            opts.Cache          = 'on';
            opts.CacheTol       = obj.tol;
            opts.TolMesh        = obj.tol;
        end
      end
      opts.MaxFunEvals    = 10000*length(x0);
      opts.MaxIter        = 10000*length(x0);
      opts.TolFun         = obj.tol;
      opts.TolX           = obj.tol;
      switch obj.alg
        case 'ga'
          problem.options = gaoptimset();
          problem.options = gaoptimset(problem.options,opts);
          problem.options = gaoptimset(problem.options,'Display',obj.display);
        case 'de'
          problem.options = opts;
          problem.options.display = obj.display;
          problem.options.lb = cons.lb;
          problem.options.ub = cons.ub;
        case 'pattern'
          problem.options = psoptimset();
          problem.options = psoptimset(problem.options,opts);
          problem.options = psoptimset(problem.options,'Display',obj.display);
        otherwise
          problem.options = opts;
          problem.options.Display = obj.display;
      end

      tic
      switch problem.solver,
        case 'fminsearch'
          [x,err,exit_flag,output] = fminsearch (eval_func, x0, problem.options);
        case 'fminunc'
          [x,err,exit_flag,output] = fminunc (eval_func, x0, problem.options);
        case 'fmincon'
          [x,err,exit_flag,output] = fmincon (eval_func, x0, problem.Aineq, problem.bineq, problem.Aeq, problem.beq, problem.lb, problem.ub, problem.nonlcon, problem.options);
        case 'ga'
          [x,err,exit_flag,output] = ga (eval_func, length(x0), problem.Aineq, problem.bineq, problem.Aeq, problem.beq, problem.lb, problem.ub, problem.nonlcon, problem.options);
        case 'de'
          [x,output] = deoptim (eval_func, length(x0), problem.options);
          exit_flag = 0;
        case 'patternsearch'
          [x,err,exit_flag,output] = patternsearch (eval_func, x0, problem.Aineq, problem.bineq, problem.Aeq, problem.beq, problem.lb, problem.ub, problem.nonlcon, problem.options);
      end
      digits = abs(round(log(obj.tol)/log(10.0))) - 1;
      x = abs(round(x*10^digits)/10^digits);
      err = eval_func (x);
      exec_time = toc;

    end
  end
end
