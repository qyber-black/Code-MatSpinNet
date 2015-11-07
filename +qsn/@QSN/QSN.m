classdef (Sealed = true) QSN < handle
  % QSN   Quantum spin network class
  %
  % This class represents a quantum spin network.
  %
  % QSN Properties:
  %   type       - type identifier
  %   N          - number of spins
  %   pos        - physical positions of spins
  %   H          - Hamiltonian
  %
  % QSN Methods:
  %   obj = QSN (type, varargin)  - C'tor
  %   eps = obj . epsilon (e)     - Get and set precision
  %   obj . disp ()               - Display quantum spin network
  %   obj . plot ()               - Plot quantum spin network
  %   p = obj . prob (t)          - Compute/get transition probabilities
  %   p = obj.trace (l, k, t)     - Calculate l-k trace
  %
  %   obj . plot_dynamics (scale, T, Pts)
  %                               - Plot spinnet dynamics
  %   [pm,tm] = obj . prob_max (tmin, tmax, steps, depth, show)
  %                               - Estimate maximum transition probability in
  %                                 time interval
  %   [p,q,theta,err,info] = obj. is_attainable (obj,ispin,jspin,x,shortest)
  %                               - Check if probability bound is attainable
  %                                 (rings only)
  %   [x,p,q,theta,OK,err,info] = find_attainable (r,ispin,jspin,s)
  %                               - Find attainability parameters
  %                                 (rings only)
  %   N = obj . opt_attainable (ispin, jspin, x, s, shortest)
  %                               - Target function for finding attainability
  %                                 (rings only; for find_attainable)
  properties (SetAccess = private)
    type   = 'undef';    % Type of quantum spin network
    N      = NaN;        % Number of spins in network
    H      = NaN;        % Network Hamiltonian
    pos    = NaN;        % Physical positions of spins
    excite = NaN;        % Number of excitations in the system
  end
  properties (Access = private)
    p_max      = NaN;                   % Maximum transition probabilities cache
    eps        = double(eps('single')); % Precision
  end
  methods
    function obj = QSN (type, varargin)
      % obj = QSN (type, varargin) - Create a quantum spin network description
      %
      %   type      - Type of spin network
      %   varargin  - Spin network decription arguments forwarded to the
      %               network type handlers:
      %               type == 'ring' - ring topology
      %                 N - Number of spins or vector with coupling strengths
      %                     between neighbouring spins
      %                 C - 'XX': XX Coupling
      %                     'H': Heisenberg coupling (optional, default: 'XX')
      %                 Z - Vector giving magnetic field strengths on
      %                     individual spins for engineered rings (optional)
      %               type == 'chain' - chain topology
      %                 N - Number of spins or vector with coupling
      %                     strengths between neighbouring spins
      %                 C - 'XX': XX Coupling
      %                     'H': Heisenberg coupling (optional, default: 'XX')
      %                 Z - Vector giving magnetic field strengths on
      %                     individual spins for engineered chains (optional)
      %               type == 'xring' - ring of chains topology
      %                 N - Number if spins in extended ring
      %                 M - Step between chains
      %                 L - Length of chains
      %                 C - 'XX': XX Coupling
      %                     'H': Heisenberg coupling (optional, default 'XX')
      %                 Z - Vector giving magnetic field strengths on individual
      %                     spins for engineered xring networks
      %               type == '2D' - 2D network with couplings from positions
      %                 N - Number if spins in 2D network
      %                 X - if string:'random', 'star', 'chain', 'xring', or '2D' topology
      %                       X == 'random' (N is total number of spins)
      %                         D - [xmin xmax; ymin max] rectangle for positions (optional, default: [0 1; 0 1])
      %                         C - 'XX': XX Coupling
      %                             'H': Heisenberg coupling (optional, default 'XX')
      %                         L - Couplings between neighbours only (optional, default 0/false)
      %                      X == 'xring' (N is number of spins in ring)
      %                         M - Step between spins in ring that have chains
      %                         L - Length of chains
      %                         C - 'XX': XX Coupling
      %                             'H': Heisenberg coupling (optional, default 'XX')
      %                         L - Couplings between neighbours only (optional, default 0/false)
      %                         E - Vector for engineered ring (optional)
      %                      X == 'chain' (N is number of spins in chain)
      %                        C - 'XX': XX Coupling
      %                            'H': Heisenberg coupling (optional, default 'XX')
      %                        L - Couplings between neighbours only (optional, default 0/false)
      %                        E - Vector for engineered chain (optional)
      %                      X == 'ring' (N is number of spins in ring)
      %                        C - 'XX': XX Coupling
      %                            'H': Heisenberg coupling (optional, default 'XX')
      %                        L - Couplings between neighbours only (optional, default 0/false)
      %                        E - Vector for engineered ring (optional)
      %                     if matrix - matrix of N given positions [x0, y0; ...; xn, yn];
      %                        C - 'XX': XX Coupling
      %                            'H': Heisenberg coupling (optional, default 'XX')
      %                        L - Couplings between neighbours only (optional, default 0/false)
      %               type == 'hes' - Higher excitation subsapce spin network
      %                 J - if 1 matrix assume isotropic Heisenberg
      %                     if 3 matrices in cell then anisotripic Heisenberg
      %                     if 2 matrices in cell assume XY interaction
      %                     if 1 matrices in cell isotropic Heisenberg
      %                 S - vector of excitation subspaces to represent,
      %                     in which order (1:k - from 1 to k <=N)
      %                 Z - vector of biases on the spins (optional)
      %   obj       - Quantum spin network object
      %
      if ~isstr(type)
        error ('Type argument must be a string');
      end
      if strcmpi (type, 'ring')
        create_uniform_1d (obj, 1, varargin);
      elseif strcmpi (type, 'chain')
        create_uniform_1d (obj, 0, varargin);
      elseif strcmpi (type, 'xring')
        create_xring (obj, varargin);
      elseif strcmpi (type, '2D')
        create_2D (obj, varargin);
      elseif strcmpi (type, 'hes')
        create_hes (obj, varargin);
      else
        error (['No constructor for spin network type ' type]);
      end
      obj.type = type;
    end
    function eps = epsilon (obj, e)
      % eps = obj . epsilon (e) - Get and set precision
      %
      %   e   - Set new precision if argument is given
      %   obj - Quantum spin network object
      %   eps - Get current precision if output variables is given
      %         (gets previous value if precision arguemnt is given)
      %
      if nargout
        eps = obj.eps;
      end
      if exist('e','var')
        if isscalar (e)
          obj.eps = e;
          obj.p_max = NaN;
        else
          warning ('Illegal precision argument ignored');
        end
      end
    end
    function disp (obj)
      % obj.disp () - Display quantum spin network
      %
      %   obj - Quantum spin network object
      %
      fprintf (1, 'Quantum Spin Network %s with %d spins, %d excitations\n\n Physical positions:\n', obj.type, obj.N, obj.excite);
      disp(obj.pos);
      fprintf (1, ' Hamiltonian:\n');
      disp(obj.H);
      fprintf (1, ' Maximum transition probabilities:');
      if (isnan(obj.p_max))
         fprintf (1, ' not computed\n\n');
      else
        fprintf (1, '\n');
        disp(obj.p_max);
      end
    end
    function plot (obj)
      % obj.plot () - Plot quantum spin network
      %
      %   obj - Quantum spin network object
      %
      clf;
      hold all;
      set(0,'defaultAxesFontName', 'Arial')
      set(0,'defaultAxesFontSize', 16)
      set(0,'defaultTextFontName', 'Arial')
      subplot(2,2,1);
      colormap(flipud(gray));
      heatmaptext (abs(obj.H), 'colorbar', false);
      title('Modulus of Hamiltonian');
      subplot(2,2,2);
      heatmaptext (angle(obj.H), 'colorbar', false, 'clim', [0 2*pi]);
      title('Angle of Hamiltonian');
      subplot(2,2,3);
      plot(obj.pos(1,:),obj.pos(2,:), 'ob');
      hold on;
      plot(obj.pos(1,:),obj.pos(2,:), '.r');
      hold off;
      xlabel('X Position');
      ylabel('Y Position');
      title('Physical positions of spins');
      subplot(2,2,4);
      obj.prob ();
      heatmaptext (obj.p_max, 'colorbar', false, 'clim', [0 1]);
      title('Maximum transition orobabilities');
      mtit(sprintf ('Quantum Spin Network %s with %d spins, %d excitations\n', obj.type, obj.N, obj.excite));
    end
    function p = prob (obj, t)
      % p = obj.prob (t) - Compute/return transition probabilities
      %
      %   t   - Time(s) for which probabilities should be computed (optional):
      %         t scalar: calculate for given time (returns matrix)
      %         t vector: calculate for given values (returns cell array of matrices)
      %         If t is not given or t is inf or nan, maximum transition probabilties
      %         are calculated / returned. This function caches the maximum transition
      %         probabilities to avoid re-calculation.
      %   obj - Quantum spin network object
      %   p   - Return value as above. Values are only returned if output variable is
      %         present. If no output value is given, maximum transition probabilities
      %         are still calculated and chaced, but other probabilty calculations are
      %         ignored.
      %
      if ~exist('t', 'var')
        t = inf;
      end
      if isscalar(t)
        if isnan(t) || isinf(t)
          compute_prob_max (obj);
          if nargout
            p = obj.p_max;
          end
        else
          if nargout
            p = abs(expm(-i*obj.H*t)).^2;
          end
        end
      elseif isvector(t)
        if nargout
          p = cell (size(t));
          for l = 1:length(t)
            p{l} = obj . prob(t(l));
          end
        end
      else
        error ('Illegal time argument');
      end
    end
    function p = trace (obj, l, k, t)
      % p = obj.trace (l, k, t) - Calculate l-k trace
      %
      %   l   - initial state
      %   k   - measured state
      %   t   - time vector
      %   obj - Quantum spin network object
      %   p   - Probabilities for l-k trace at times t
      %
      [V,e] = eig (obj.H);
      e = diag(e);
      E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (t), 'UniformOutput', false);
      p =  cellfun (@(x) abs(x(k, l))^2, E);
    end
  end
end
