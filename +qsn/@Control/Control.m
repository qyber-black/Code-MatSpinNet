classdef (Sealed = true) Control < handle
  % Control   Quantum spin network control problem class
  %
  % This class represents a quantum spin network control problem.
  %
  % Control Properties:
  %   system               - QSN object
  %   problem              - Control problem description
  %   H                    - Cell array of Hamiltonians
  %
  % Control Methods:
  %   obj = Control (s, type, varargin)   - C'tor
  %   obj . disp ()                       - Display quantum spin network control problem
  %   [x,err,exec_time,exit_flag,output,x0] = obj.solve_switch (opt, x0)
  %                                       - solve problem by switching controls on and off
  %   [x,err,exec_time,exit_flag,output,x0] = obj.solve_static (opt, x0)
  %                                       - solve problem by static controls
  %

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  properties (SetAccess = private)
    system    = 'undef';                   % QSN object
    problem   = struct ('type', 'undef');  % Control problem description
    H         = NaN;                       % Cell array of Hamiltonians
  end
  properties (Access = private)
    H_N        = NaN;                       % Number of controls
    E          = NaN;                       % Eigenvalues of Hamiltonians
    V          = NaN;                       % Eigenvectors of Hamiltonians
    Vab        = NaN;                       % Products of consecutive eigenvector matrices
  end
  methods
    function obj = Control (s, type, ctrls, varargin)
      % obj = Control (s, type, varargin) - C'tor
      %
      %   s     - QSN object to be controlled
      %   ctrls - Cell array of control Hamiltonians
      %   type  - Type of control problem
      %           type == 'state_transfer'
      %             varargin{1} : init     - initial state
      %             varargin{2} : target   - target state
      %             varargin{3} : time     - target time (set to NaN for undef, static only)
      %             varargin{4} : max_bias - maximum value for bias (set NaN for undef, static only)
      %             varargin{5} : dt       - readout time (inverse bandwidth in units of J) [optional]
      %             varargin{6} : pertubation - pertubation of drift hamiltonian
      %             varargin{7} : samples     - samples for averaging over perturbation
      %           type == 'operator'
      %             varargin{1} : u_taret - target operator
      %             varargin{2} : time    - target time
      %
      % FIXME: helper constructs to setup ctrls?
      if ~isa (s, 'qsn.QSN')
        error ('Illegal system to control');
      end
      obj.system = s;
      if iscell (ctrls)
        obj.H = cell (1,length (ctrls) + 1);
        obj.H{1} = obj.system.H;
        for l = 1:length(ctrls)
          if ~ismatrix(ctrls{l}) || size(ctrls{l},1) ~= obj.system.N || size(ctrls{l},2) ~= obj.system.N
            error ('Illegal control Hamiltonian');
          else
            obj.H{l+1} = ctrls{l};
          end
        end
      elseif ismatrix (ctrls)
        if size(ctrls,1) ~= obj.system.N || size(ctrls,2) ~= obj.system.N
          error ('Illegal control Hamiltonian');
        end
        obj.H = cell(1,2);
        obj.H{1} = obj.system.H;
        obj.H{2} = ctrls;
      else
        error ('Illegal control Hamiltonians');
      end
      if strcmpi (type, 'state_transfer')
        obj.problem = struct ('type',     'state_transfer', ...
                              'init',     varargin{1}, ...
                              'target',   varargin{2}, ...
                              'time',     varargin{3}, ...
                              'max_bias', varargin{4});
        if size(varargin,2) > 4
          obj.problem.dt = varargin{5};
        end
        if size(varargin,2) > 5
          obj.problem.perturbation = varargin{6};
        end
        if size(varargin,2) > 6
          obj.problem.samples = varargin{7};
        end
        if obj.problem.init ~= double(uint64(obj.problem.init)) || obj.problem.init < 1 || obj.problem.init > obj.system.N
          error ('Illegal initial state');
        end
        if obj.problem.target ~= double(uint64(obj.problem.target)) || obj.problem.target < 1 || obj.problem.target > obj.system.N
          error ('Illegal target state');
        end
        if ~isscalar(obj.problem.time) || obj.problem.time < 0.0
          error ('Ilelgal target time');
        end
      elseif strcmpi (type, 'operator')
        obj.problem = struct ('type',    'operator', ...
                              'operator', varargin{1}, ...
                              'time',   varargin{2})
        if ~ismatrix(obj.problem.operator) || size(obj.problem.operator, 1) ~= obj.system.N || size(obj.problem.operator,2) ~= obj.system.N
          error ('Illegal target operator');
        end
        if ~isscalar(obj.problem.time) || obj.problem.time < 0.0
          error ('Ilelgal target time');
        end
      else
        error ('Unknown control problem type');
      end
    end
    function [err,grad,prob] = test_static_st_dt_avg (obj, x)
      % [err,grad,prob] = obj.test_static_st_avg (x)
      %
      % Test target function for static state transfer problems
      % with time window, average over pertubations using sampling
      %
      %   x       - Value
      %   obj     - Control object
      %   err     - Error
      %   prob    - Probability error
      %
      if isnan(obj.problem.time)
        T = abs(x(end)); % ensure T >=0
      else
        T = obj.problem.time;
      end
      H = obj.H{1};
      for l = 2:obj.H_N
        H = H + x(1,l-1) * obj.H{l};
      end
      N = size(H,1);
      out = zeros(N,1);
      in  = zeros(N,1);
      out(obj.problem.target)= 1;
      in(obj.problem.init)   = 1;
      prob = 0;
      ppp=[];
      for rep = 1:obj.problem.samples*10
        HP = obj.H{1} .* (obj.problem.perturbation * (rand(size(H)) * 2 - 1));
        HP = triu(HP) + triu(HP,1)';
        [V,E] = eig(H + HP);
        E = diag(E);
        N = length(E);
        L = ones(N,1)*E.'-E*ones(1,N);
        W = (out'*V).*(V'*in)';
        pp = 1-real(W*(exp(1i*L*T).*sinc(L*obj.problem.dt/(2*pi)))*W');
        prob = prob + pp;
        ppp(rep) = pp;
      end
      prob = prob / obj.problem.samples;
      err = prob / obj.problem.samples;
      grad = [];
      NN = 1:obj.problem.samples*10;
      scatter(NN,ppp,3,ppp,'o','filled');
      hold on;
      plot(NN,cumsum(ppp)./NN,'b');
      ylabel('Probability error');
      xlabel('Samples');
      title('Probability error caused by Hamiltonian pertubation using sampling and average');
      hold off;
    end
    function disp (obj)
      % obj . disp () - Display quantum spin network control problem
      %
      %   obj - Quantum spin network control problem object
      %
      fprintf (1, 'System to control:\n\n');
      disp(obj.system);
      if strcmp (obj.problem.type, 'state_transfer')
        fprintf (1, 'State transfer control problem from node %d to %d at time %g\n\n', obj.problem.init, obj.problem.target, obj.problem.time);
      elseif strcmp (obj.problem.type, 'operator')
        fprintf (1, 'Target operator control problem at time %g: \n\n', obj.problem.time);
        disp(obj.problem.operator);
      else
        error ('Unknown control problem type');
      end
      fprintf (1, 'With %d Hamiltonians: ', length(obj.H));
      disp (obj.H);
    end
    function [x,err,exec_time,exit_flag,output,x0] = solve_switch (obj, opt, x0)
      % [x,err,exec_time,exit_flag,output,x0] = obj.solve_switch (opt, x0) - solve problem by switching controls on and off
      %
      % System Hamiltonian alone or system Hamiltonian and one of the control Hamiltonians can be on.
      %
      %   opt       - Opt object used for optimisation
      %   x0        - Number of intial switch times selected randomly or
      %               vector of initial switching times
      %   obj       - Control object
      %   x         - Switching times (x(n) is time Hamiltonian x%H_N is on, incl. system Hamiltonian 0)
      %   err       - Error
      %   exec_time - Execution time
      %   exit_flag - Optimisation exit flag status indicator
      %   output    - Info about optimisation
      %   x0        - Initial switching times
      %
      % FIXME: report solution suitably
      obj.H_N = length(obj.H);
      % Initial switch times
      if length (x0) == 1
        x = rand (x0, obj.H_N);
        x = x ./ (ones(x0,1) * sum (x, 1)) * obj.problem.time;
        x0 = reshape (x', x0 * obj.H_N, 1);
      end
      switches = fix(length (x0) / obj.H_N);
      % Add system Hamiltonian to control Hamiltonians
      HH{1} = obj.H{1};
      HH{2:obj.H_N} = obj.H{2:end} + obj.H{1};
      % Compute eigendecomposition of all Hamiltonians
      [obj.V,e] = cellfun(@(x) eig(x), HH, 'UniformOutput', false);
      obj.E = cellfun(@(x) diag(x), e, 'UniformOutput', false);
      % Define objective and constraints
      if strcmpi (obj.problem.type, 'state_transfer')
        eval_func = @(x,target,digits) obj.eval_sw_st (x);
      elseif strcmpi (obj.problem.type, 'operator')
        eval_func = @(x,target,digits) obj.eval_sw_op (x);
        % Compute Vm'*V_{m+1}
        obj.Vab = cellfun(@(x,y) x'*y, obj.V, {obj.V{2:end}, obj.V{1}}, 'UniformOutput', false);
      else
        error ('Unknown control problem');
      end
      if ~isa (opt, 'qsn.Opt')
        error ('Illegal optimiser');
      end
      [x,err,exec_time,exit_flag,output] = opt.solve (eval_func, x0);
    end
    function [x,err,exec_time,exit_flag,output,x0] = solve_static (obj, opt, x0)
      % [x,err,exec_time,exit_flag,output,x0] = obj.solve_static (opt, x0) - solve problem via static control
      %
      %   opt       - Opt object used for optimisation
      %   x0        - Initial values for switches (if empty, selected randomly)
      %   obj       - Control object
      %   x         - Optimal biases
      %   err       - Error
      %   exec_time - Execution time
      %   exit_flag - Optimisation exit flag status indicator
      %   output    - Info about optimisation
      %   x0        - Initial biases
      %
      % FIXME: report solution suitably
      obj.H_N = length(obj.H);
      % Initial biases
      p1 = obj.problem.init;
      p2 = obj.problem.target;
      T = obj.problem.time;
      if isempty(T) || isnan(T)
        obj.problem.time = NaN;
        T = rand(1)*20*(obj.H_N-1);
      end
      if ~exist('x0', 'var') || isempty(x0)
        x0 = rand (1,obj.H_N-1);
        if isnan(obj.problem.time)
          x0 = [x0 T];
        end
      end

      % Define objective and constraints
      if strcmpi (obj.problem.type, 'state_transfer')
        if isfield(obj.problem,'dt')
          if isfield(obj.problem,'perturbation')
            eval_func = @(x) obj.eval_static_st_dt_avg(x);
          else
            eval_func = @(x) obj.eval_static_st_dt(x);
          end
        else
          eval_func = @(x) obj.eval_static_st(x);
        end
      else
        error ('Unknown control problem');
      end
      if ~isa (opt, 'qsn.Opt')
        error ('Illegal optimiser');
      end

      if strcmp(opt.alg,'de')
        % DE optimisation
        if isnan(obj.problem.max_bias)
          obj.problem.max_bias = 1000;
        end
        cons.Aineq     = [];
        cons.bineq     = [];
        cons.Aeq       = [];
        cons.beq       = [];
        cons.nonlcon   = [];
        cons.lb = zeros(1,length(x0));
        cons.ub = obj.problem.max_bias * ones(1,length(x0));
        if isnan(obj.problem.time)
          cons.ub(1,length(x0)) = 20 * (obj.H_N-1);
        end
        [x,err,exec_time,exit_flag,output] = opt.solve (eval_func, x0, cons);
      else
        % Other Optimisation
        if ~isnan(obj.problem.max_bias)
          lx = length(x0);
          cons.Aineq     = [];
          cons.bineq     = [];
          cons.Aeq       = [];
          cons.beq       = [];
          % cons.lb        = zeros(lx,1);
          cons.lb        = -obj.problem.max_bias*ones(lx,1);
          cons.ub        = obj.problem.max_bias*ones(lx,1);
          if isnan(obj.problem.time)
            cons.ub(lx,1) = 100;
          end
          cons.nonlcon   = [];
          [x,err,exec_time,exit_flag,output] = opt.solve (eval_func, x0, cons);
        else
          [x,err,exec_time,exit_flag,output] = opt.solve (eval_func, x0);
        end
      end
      [~,~,err] = eval_func(x); % Ensure fidelity error is reported
    end
    % FIXME: other solve functions
  end
  methods %(Access=private)
    function [err,grad] = eval_sw_op (obj, x)
      % [err,grad] = obj . eval_sw_op (x)
      %
      % Target function for operator problems
      %
      %   x       - Value
      %   obj     - Control object
      %   err     - Error
      %   grad    - Gradient
      %
      % FIXME: not verified
      x = abs(x);
      tm = cellfun(@(n) x(n:obj.H_N:end), num2cell([1:obj.H_N]), 'UniformOutput', false);
      for m = 1:obj.H_N
        D{m} = cellfun (@(x) diag(exp(-i*x*obj.E{m})), num2cell(tm{m}), 'UniformOutput', false);
        U{m} = cellfun (@(x) obj.V{m}*x*obj.V{m}', D{m}, 'UniformOutput', false);
      end
      N  = length (obj.E{1});
      Uf = obj.V{1};
      for k=1:fix(length(x)/obj.H_N)
        for m=1:obj.H_N
          Uf = Uf * D{m}{k} * obj.Vab{m};
        end
      end
      Uf = Uf * obj.V{1}';
      err = 0.5 * (1 - real(trace(obj.problem.opreator' * Uf)) / N);
      if nargout > 1,
        for m=1:obj.H_N
          U{m} = cellfun(@(x) obj.V{m} * x * obj.V{m}', D{m}, 'UniformOutput', false);
        end
        U0 = eye (N);
        for k=1:(length(x)/obj.H_N)
          for m=1:obj.H_N
            U0 = U0 * U{m}{k};
            Uf = U{m}{k}' * Uf;
            dU = -i * U0 * obj.H{m} * Uf;
            grad(m,k) = -0.5*imag(trace(obj.problem.operator' * dU))/N;
          end
        end
        grad = reshape(grad,length(x),1);
      end
    end
    function [err,grad,prob] = eval_static_st (obj, x)
      % [err,grad,prob] = obj.eval_static_st (x)
      %
      % Target function for static state transfer problems
      %
      %   x       - Value
      %   obj     - Control object
      %   err     - Error
      %   prob    - Probability error
      %
      if isnan(obj.problem.time)
        T = abs(x(end)); % ensure T >=0
      else
        T = obj.problem.time;
      end
      H = obj.H{1};
      for l = 2:obj.H_N
        H = H + x(1,l-1) * obj.H{l};
      end
      %U = expm(-1i*T*H);
      [V,D] = eig(-1i*T*H);
      U = V * diag(exp(diag(D))) / V;
      %
      phi = U(obj.problem.target,obj.problem.init);
      prob = 1-abs(phi)^2; % FIXME: minimise time?
      err = prob;
      if nargout > 1,
        grad = zeros (obj.H_N-1,1); % derivative of the error wrt f(mk)
        for m=1:obj.H_N-1
          [~,dU]  = dexpma(-1i*T*H,-1i*T*obj.H{m+1});
          grad(m) = -2*real(dU(obj.problem.target,obj.problem.init)*conj(phi));
        end
        % dF/dT = -2 Re <p2|-iH*U|p1> conj(<p2|U|p1>)
        if isnan(obj.problem.time)
          grad(obj.H_N) = -2*imag(H(obj.problem.target,:)*U(:,obj.problem.init)*conj(phi));
        end
      end
    end
    function [err,grad,prob] = eval_static_st_dt (obj, x)
      % [err,grad,prob] = obj.eval_static_st_dt (x)
      %
      % Target function for static state transfer problems
      % with time window
      %
      %   x       - Value
      %   obj     - Control object
      %   err     - Error
      %   prob    - Probability error
      %
      if isnan(obj.problem.time)
        T = abs(x(end)); % ensure T >=0
      else
        T = obj.problem.time;
      end
      H = obj.H{1};
      for l = 2:obj.H_N
        H = H + x(1,l-1) * obj.H{l};
      end
      N = size(H,1);
      out = zeros(N,1);
      in  = zeros(N,1);
      out(obj.problem.target)= 1;
      in(obj.problem.init)   = 1;
      [V,E] = eig(H);
      E = diag(E);
      N = length(E);
      L = ones(N,1)*E.'-E*ones(1,N);
      W = (out'*V).*(V'*in)';
      prob = 1-real(W*(exp(1i*L*T).*sinc(L*obj.problem.dt/(2*pi)))*W');
      err = prob;
      grad = [];
    end
    function [err,grad,prob] = eval_static_st_dt_avg (obj, x)
      % [err,grad,prob] = obj.eval_static_st_avg (x)
      %
      % Target function for static state transfer problems
      % with time window, average over pertubations using sampling
      %
      %   x       - Value
      %   obj     - Control object
      %   err     - Error
      %   prob    - Probability error
      %
      if isnan(obj.problem.time)
        T = abs(x(end)); % ensure T >=0
      else
        T = obj.problem.time;
      end
      H = obj.H{1};
      for l = 2:obj.H_N
        H = H + x(1,l-1) * obj.H{l};
      end
      N = size(H,1);
      out = zeros(N,1);
      in  = zeros(N,1);
      out(obj.problem.target)= 1;
      in(obj.problem.init)   = 1;
      prob = 0;
      for rep = 1:obj.problem.samples
        HP = obj.H{1} .* (obj.problem.perturbation * (rand(size(H)) * 2 - 1));
        HP = triu(HP) + triu(HP,1)';
        [V,E] = eig(H + HP);
        E = diag(E);
        N = length(E);
        L = ones(N,1)*E.'-E*ones(1,N);
        W = (out'*V).*(V'*in)';
        pp = 1-real(W*(exp(1i*L*T).*sinc(L*obj.problem.dt/(2*pi)))*W');
        prob = prob + pp;
      end
      prob = prob / obj.problem.samples;
      err = prob / obj.problem.samples;
      grad = [];
    end
    function [err,grad] = eval_sw_st (obj, x)
      % [err,grad] = obj.eval_sw_st (x)
      %
      % Target function for state transfer problems
      %
      %   x       - Value
      %   obj     - Control object
      %   err     - Error
      %   grad    - Gradient
      %
      x = abs(x);
      U = cell(obj.H_N,1);
      parfor m = 1:obj.H_N
        U{m} = cellfun(@(v) obj.V{m}*diag(exp(-i*v*obj.E{m}))*obj.V{m}', num2cell(x(m:obj.H_N:end)), 'UniformOutput', false);
      end
      N  = length (obj.E{1});
      Uf = eye (N);
      for k = 1:fix(length(x)/obj.H_N)
        for m = 1:obj.H_N
          Uf = Uf*U{m}{k};
        end
      end
      C = Uf(obj.problem.target, obj.problem.init);
      err = 1 - abs (C)^2;
      if nargout > 1,
        U0 = eye(N);
        for k = 1:fix (length(x)/obj.H_N)
          for m = 1:obj.H_N
            U0 = U0 * U{m}{k};
            Uf = U{m}{k}' * Uf;
            dU = U0 * obj.H{m} * Uf;
            grad(m,k) = -2*(imag(dU(obj.problem.target,obj.problem.init)) * real(C) - real(dU(obj.problem.target,obj.problem.init))*imag(C));
          end
        end
        grad = reshape (grad,length(x),1);
      end
    end
  end
end
