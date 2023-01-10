classdef (Sealed = true) DS < handle
  % DS   Distance space class
  %
  % This class represents a discrete distance space.
  %
  % DS Properties:
  %   N     - number of nodes
  %   D     - Matrix of distances between nodes
  %
  % DS Methods:
  %   obj = DSN (d, varargin)              - Create a distance space
  %   eps = obj . epsilon (e)              - Get and set precision
  %   obj . disp ()                        - Display distance space
  %   obj . plot ()                        - Plot distance space
  %   I = obj . inertia (alpha)            - Compute inertia
  %   [delta4,scaled_delta4,g] = obj . Gromov4pt(ind)
  %                                        - Compute Gromove 4-point deltas
  %   [max_delta4,max_scaled_delta4] = obj . max_Gromov4pt () 
  %                                        - Calculate maximum Gromov 4-point deltas
  %   [x,I] = obj . check_non_negative ()  - Check for negative distances
  %   [x,I] = obj . check_coincidence ()   - Check for coincident nodes
  %   [x,I] = obj . check_symmetry ()      - Check for symmetric distances
  %   [x,Tmin,T] = obj . check_triangle_inequality (v)
  %                                        - Check if distances fulfill triangle inequality
  %   [x,n,c,s,t] = obj . check_metric (v) - Check if distances form a metric
  %   res = obj . generate_metric (symmetry, coincidence, triangle, tmod, verbose)
  %                                        - Fix metric violations and generate a new distance
  %                                          space object
  %   diam = obj . diameter ()             - Diameter of distance space
  %   lambda = obj . gramsaddle (kappa)    - Are distances embeddable in negative curvature space?
  %   lambda = obj . gramsphere (R)        - Are distances embeddable in a sphere?
  %   lambda = obj . cmflat ()             - Are distances embeddable in Euclidean space?
  %   [G,lambda]  = obj . gram_matrix ()   - Calculate Gram matrix
  %
  %   obj . plot_graph (type, scale, violations) 
  %                                        - Plot distances as graph

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  properties (SetAccess = private)
    N    = NaN;        % Number of nodes
    D    = NaN;        % Distance matrix
  end
  properties (Access = private)
    diam = NaN;                   % Cache for distance space diameter
    eps  = double(eps('single')); % Precision
  end
  methods
    function obj = DS (d, varargin)
      % obj = DSN (d, varargin) - Create a distance space
      %
      %   d        - Matrix of distances between nodes or
      %              qsn.QSN object
      %   varargin - If d is qsn.QSN object:
      %              {1}: t-     Time for which probabilities are used to generate distance sapce
      %                          (optional, default: inf)
      %              {2}: func - Functional handle argument to compute distance from probability
      %                          (optional, default: @(x) -log(x))
      %   obj  - Distance space object
      %
      if isa (d, 'qsn.QSN')
        if size(varargin,2) == 0
          t = inf;
        else
          t = varargin{1};
          if ~isscalar (t)
            error ('Illegal time value');
          end
        end
        if size(varargin,2) < 2
          func = @(x) -log(x);
        else
          func = varargin{2};
          if ~isa(func, 'function_handle')
            error ('Illegal function handle');
          end
        end
        obj.D = func (d.prob (t));
        obj.N = size (obj.D, 1);
      else
        if ~ismatrix(d) || diff(size(d)) ~= 0
          error ('Distances must be given as a square matrix');
        end
        obj.D = d;
        obj.N = size (d, 1);
      end
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
        else
          warning ('Illegal precision argument ignored');
        end
      end
    end
    function disp (obj)
      % obj . disp () - Display distance space
      %
      %   obj - Distance space object
      %
      fprintf (1, 'Distance space with %d nodes\n\n Distances:\n', obj.N);
      disp (obj.D);
    end
    function plot (obj)
      % obj . plot () - Plot distance space
      %
      %   obj - Distance space object
      %
      clf;
      hold all;
      set(0,'defaultAxesFontName', 'Arial')
      set(0,'defaultAxesFontSize', 16)
      set(0,'defaultTextFontName', 'Arial')
      subplot(2,2,1);
      colormap(flipud(gray));
      heatmaptext (round(obj.D / obj.eps)*obj.eps, 'colorbar', false);
      freezeColors ();
      title('Distances');
      subplot(2,2,2);
      colormap('default');
      cmap = colormap ();
      [mx my] = size(cmap);
      I = obj.inertia();
      bar_h = bar(1:obj.N,I);
      bar_ch = get (bar_h, 'Children');
      set (bar_ch, 'EdgeColor', 'none');
      set (bar_ch, 'CDataMapping', 'direct');
      min_I = min(I);
      max_I = max(I);
      if max_I - min_I < obj.eps
        ci = floor(mx/2);
      else
        ci = floor((I - min_I) / (max_I - min_I) * (mx-1)) + 1;
      end
      set (bar_ch, 'CData', ci);
      freezeColors ();
      axis([0 obj.N+1 0 max_I*1.05]);
      xlabel('Spin Index');
      ylabel('Inertia');
      title('Inertia, $\alpha = 2$', 'interpreter','latex','fontsize',16);
      subplot(2,2,3);
      [delta4,scaled_delta4] = obj.Gromov4pt();
      plot(delta4);
      ms = max(delta4);
      hold on;
      plot ([0 length(delta4)], [ms ms], 'r-');
      text (0.05*length(delta4), ms*1.02, sprintf('Maximum %g',ms));
      axis ([0 length(delta4) 0 ms*1.1]);
      hold off;
      xlabel('Quadrilateral Index');
      ylabel('$\delta_4$','interpreter','latex','fontsize',16);
      title('Gromov $\delta_{4}$','interpreter','latex','fontsize',16);
      subplot(2,2,4);
      plot(scaled_delta4);
      ms = max(scaled_delta4);
      hold on;
      plot ([0 length(scaled_delta4)], [ms ms], 'r-');
      text (0.05*length(scaled_delta4),ms*1.02,sprintf('Maximum %g',ms));
      axis ([0 length(scaled_delta4) 0 ms*1.1]);
      hold off;
      xlabel('Quadrilateral Index');
      ylabel('Scaled $\delta_4$','interpreter','latex','fontsize',16);
      title('Scaled Gromov $\delta_{4}$','interpreter','latex','fontsize',16)
      mtit(sprintf ('Distance space with %d nodes\n', obj.N));
    end
    function I = inertia (obj, alpha)
      % I = obj . inertia (alpha) - Compute inertia \sum_i |d(i,j)|^alpha
      %
      %   alpha - Exponent (optional, default: 2.0)
      %   obj   - Distance space object
      %   I     - Vector of intertia values of nodes
      %
      if ~exist('alpha','var')
        alpha = 2.0;
      end
      I = sum(abs(obj.D).^alpha);
    end
    function [delta4,scaled_delta4,G] = Gromov4pt(obj)
      % [delta4,scaled_delta4,g] = obj . Gromov4pt () - Calculate Gromove 4-point deltas
      %
      % 4-point Gromov \delta (delta1) and scaled 4-point gromove \delta (delta2) for all
      % all quadrangles in the distance space. G contains the sorted "diameter sums" of all
      % quadrangles.
      %
      % With sums of lenghts of opposite diagonals
      %    L >= M >= S
      % For infinite graph:
      %    \delta_G(G) = sup(L(QUAD) - M(QUAD)) / 2 : QUAD \in G} < \infty
      % For finite graph scale by diameter:
      %    \delta_G(QUAD) / diam(QUAD) < (\sqrt(2) - 1) / \sqrt(2) ~ .29
      %
      %
      %   obj           - Distance space object
      %   delta4        - delta_4 vector
      %   scaled_delta4 - scaled delta_4 vector
      %   g             - sorted "diameter sums" of all quadrangles
      %
      if obj.N < 4
        error ('Gromov4pt requires at least four nodes');
      end
      ind = nchoosek([1:obj.N],4);
      S = [obj.N obj.N];
      ind12 = sub2ind(S,ind(:,1),ind(:,2));
      ind34 = sub2ind(S,ind(:,3),ind(:,4));
      ind13 = sub2ind(S,ind(:,1),ind(:,3));
      ind24 = sub2ind(S,ind(:,2),ind(:,4));
      ind14 = sub2ind(S,ind(:,1),ind(:,4));
      ind23 = sub2ind(S,ind(:,2),ind(:,3));
      G = sort ([obj.D(ind12) + obj.D(ind34), obj.D(ind13) + obj.D(ind24), obj.D(ind14) + obj.D(ind23)], 2);
      delta4 = (G(:,3)-G(:,2))/2;
      scaled_delta4 = delta4./sum(G,2);
    end
    function [max_delta4,max_scaled_delta4] = max_Gromov4pt(obj)
      % [max_delta4,max_scaled_delta4] = obj . max_Gromov4pt () - Calculate maximum Gromov 4-point deltas
      %
      %   obj           - Distance space object
      %   delta4        - delta_4 vector
      %   scaled_delta4 - scaled delta_4 vector
      %
      [d4,sd4] = obj.Gromov4pt();
      max_delta4 = max(d4);
      max_scaled_delta4 = max(sd4);
    end
    function [x,I] = check_non_negative (obj)
      % [x,I] = obj . check_non_negative () - Check for negative distances
      %
      %   obj - Distance space object
      %   x   - Boolean indicating if all distances are non-negative
      %         (within precision)
      %   I   - Matrix mask giving non-negative distances
      %
      I = obj.D < -obj.eps;
      x = sum(sum(I)) == 0;
    end
    function [x,I] = check_coincidence (obj)
      % [x,I] = obj . check_coincidence () - Check for coincident nodes
      %
      %   obj - Distance space object
      %   x   - Boolean indicating if there are no coincident nodes
      %         and equal nodes are coincident
      %         (within precision)
      %   I   - Matrix mask giving coincident nodes
      %         (off-diagonal zero distances and diagonal zero distances)
      %
      I = (abs(obj.D) < obj.eps) .* (ones(obj.N,obj.N) - eye (obj.N));
      I = I + diag(abs(diag(obj.D)) >= obj.eps);
      x = sum(sum(I)) == 0;
    end
    function [x,I] = check_symmetry (obj)
      % [x,I] = obj . check_symmetry () - Check for symmetric distances
      %
      %   obj - Distance space object
      %   x   - Boolean indicating if the distances are symmetric
      %         (within precision)
      %   I   - Matrix mask giving non-symmetric nodes
      %
      I = abs (obj.D - obj.D') > obj.eps;
      x = sum(sum(I)) == 0;
    end
    function [x,Tmin,T] = check_triangle_inequality (obj, verbose)
      % [x,Tmin,T] = obj . check_triangle_inequality (v) - Check if distances fulfill triangle inequality
      %
      %   v     - Verbose output: plot triangle inequality differences / violations
      %           (optional, default: false)
      %   obj   - Distnace space object
      %   x     - Boolean indicating if triangle inequality is fulfilled
      %   Tmin  - Smallest triangle inequality difference
      %   T     - Triangle inequality differences (3D array)
      T = zeros(obj.N,obj.N,obj.N);
      for k=1:obj.N
        for l=1:obj.N
          for j=1:obj.N
            T(k,l,j) = obj.D(k,l) + obj.D(l,j) - obj.D(j,k);
          end
        end
      end
      Tmin = min(min(min(T)));
      err = sum(sum(sum(T<-obj.eps)));
      x = err < obj.eps;
      if exist('verbose', 'var') && verbose
        clf;
        hold all;
        set(0,'defaultAxesFontName', 'Arial')
        set(0,'defaultAxesFontSize', 28)
        set(0,'defaultTextFontName', 'Arial')
        err = sort(-reshape(T.*(T<-obj.eps), prod(size(T)), 1), 1, 'descend');
        plot ([1:length(err)],sort (err, 1, 'descend'), 'LineWidth', 2);
        title ('Triangle violations');
        xlabel ('Triangles sorted by violation');
        ylabel ('Triangle violation');
        hold off;
      end
    end
    function [x,n,c,s,t] = check_metric (obj, v)
      % [x,n,c,s,t] = obj . check_metric (v) - Check if distances form a metric
      %
      %   v   - If true, display test results (optional, default: false)
      %   obj - Distance space object
      %   x   - Boolean indicating if distnaces are a metric
      %   n   - Boolean indicating if all distances are non-negative
      %   c   - Boolean indicating if there are no coincident nodes
      %         and equal nodes are coincident
      %   s   - Boolean indicating if the distances are symmetric
      %         (within precision)
      %   t   - Boolean indicating if triangle inequality is fulfilled
      %
      n = obj.check_non_negative ();
      c = obj.check_coincidence ();
      s = obj.check_symmetry ();
      t = obj.check_triangle_inequality ();
      x = n && c && s && t;
      if exist ('v', 'var') && v
        if n
          display ('Non-negativity:      OK');
        else
          display ('Non-negativity:      violated');
        end
        if c
          display ('Coincidence:         OK');
        else
          display ('Coincidence:         violated');
        end
        if s
          display ('Symmetry:            OK');
        else
          display ('Symmetry:            violated');
        end
        if t
          display ('Triangle inequality: OK');
        else
          display ('Triangle inequality: violated');
        end
      end
    end
    function res = generate_metric (obj,symmetry,coincidence,triangle,tmode,verbose)
      % res = obj . generate_metric (symmetry, coincidence, triangle, tmod, verbose)
      %
      % Fix metric violations and generate a new distance space object.
      %
      % Create metric space from distance space
      %   symmetry    - How to fix symmetry: 'min', 'max', 'mean'
      %   coincidence - How to fix coindicence: 'min', 'max', 'mean'
      %   triagnle    - How to fix triangle inequality: 'isoceles', 'straight', 'mincon', 'minuncon'
      %   tmode       - Triangle fix mode: @(x)mean(x), @(x)min(x), @(x)max(x) (only used for isoceles)
      %   verbose     - Boolean indicating whether to report progress (optional, default: false)
      %   res         - New distance space object
      %
      if ~obj.check_non_negative () % Non-negativity
        error ('Cannot fix metric space with negative distances');
      end
      % Round distances to precision
      D = round (obj.D / obj.eps) * obj.eps;
      % ...and set (small) negative distances to zero.
      D(find(D < 0)) = 0;
      % Fix symmetry
      if ~obj.check_symmetry ()
        if verbose
          fprintf (1, 'Fixing symmetry violations\n');
        end
        if strcmpi(symmetry, 'mean')
          D = (D + D') / 2;
        elseif strcmpi(symmetry, 'min')
          D = min(D, D');
        elseif strcmpi(symmetry, 'max')
          D = max(D, D');
        else
          error ('Unknown symmetry resolution technique');
        end
      end
      % Coincidence identification
      [x,I] = obj.check_coincidence ();
      if ~x
        % Identify points
        if verbose
          fprintf (1, 'Fixing coincidence violations\n');
        end
        if I ~= I'
          error ('Cannot fix distances: coincidences are not symmetric');
        end
        [R,C] = find(I);
        if strcmpi(coincidence,'mean')
          D(R,:) = (D(R,:) + D(C,:)) / 2;
        elseif strcmpi(coincidence,'min')
          D(R,:) = min(D(R,:),D(C,:));
        elseif strcmpi(coincidence,'max')
          D(R,:) = max(D(R,:),D(C,:));
        else
          error ('Unknown coincidence resolution');
        end
        D(:,C) = D(R,:)';
        id = [1:obj.N];
        id(C(find(C > R))) = 0;
        id = find(id > 0);
        D = D(id,id);
        % ...and ensure points are coincident to themsleves
        D = D .* (ones(size(D)) - eye(size(D)));
      end
      % Fix triangle inequality
      if ~obj.check_triangle_inequality ()
        if verbose
          fprintf (1, 'Fixing triangle inequality violations\n');
        end
        if strcmpi(triangle, 'isoceles')
          D = fix_triangle_violations_isoceles (D, obj.eps, tmode, verbose);
        elseif strcmpi(triangle, 'straight')
          D = fix_triangle_violations_straight (D, obj.eps, verbose);
        elseif strcmpi(triangle, 'mincon')
          D = fix_triangle_violations_mincon (D, 2, 0, 2, verbose);
        elseif strcmpi(triangle, 'minuncon')
          D = fix_triangle_violations_minuncon (D, obj.eps, 2, verbose);
        else
          error ('Unknown triangle inequlaity correction method');
        end
      end
      % Return result
      res = qsn.DS (D);
      res.epsilon (obj.eps);
    end
    function diam = diameter(obj)
      % diam = obj . diameter () - Diameter of distance space
      %
      %   obj  - Distance space object
      %   diam - Diameter of distance space (maximum distance, chached)
      %
      if (strcmp(obj.diam,'undef'))
        obj.diam = max(max(obj.D));
      end
      diam = obj.diam;
    end
    function lambda = gramsaddle (obj,kappa)
      % lambda = obj . gramsaddle (kappa) - Are distances embeddable in negative curvature space?
      %
      % This function checks whether the distance data D is embeddable in a negatively curved
      % space of curvature kappa < 0. All eigenvalues lambda of Gram matrix
      % cosh(D*sqrt(-Kappa)) must be negative except one that has to be positive.
      %
      %   kappa  - curvature of space (positive value for negative curvature)
      %   obj    - Distance space object
      %   lambda - Eigenvalues of Gram matrix
      %
      lambda = sort(eig(cosh(D*sqrt(-kappa))));
    end
    function lambda = gramsphere (obj,R)
      % lambda = obj . gramsphere (R) - Are distances embeddable in a sphere?
      %
      % This function checks whether the distance data D is embeddable in a
      % sphere of radius R. All eigenvalues lambda of Gram matrix
      % cos(D*sqrt{Kappa}) must be positive.
      %
      %   R      - Radius of sphere
      %   obj    - Distance space object
      %   lambda - Eigenvalues of Gram matrix
      %
      if obj.diameter () > R * pi
        error ('diameter condition violated');  % diam <= pi/sqrt{K} is violated
      end
      lambda = sort(eig(cos(obj.D / R)));
    end
    function lambda = cmflat(obj)
      % lambda = obj . cmflat () - Are distances embeddable in Euclidean space?
      %
      % This function check embeddability of the distance data D in Euclidean
      % space by resorting to the Cayley-Menger matrix. The function returns
      % the eigenvalue in increasing order. All eigenvalues of Cayley-Menger
      % matrix must be negative, except one that must be positive.
      %
      %   obj    - Distance spcae object
      %   lambda - Eigenvalues of the Cayley-Menger matrix
      %

      % Cayley-Menger matrix
      CM = [0 ones(1,obj.N); ones(obj.N,1) (obj.D).^2]
      % Check eigenvalues
      lambda = sort(real(eig(CM)));
    end
    function [G,lambda] = gram_matrix (obj)
      % [G,lambda] = obj . gram_matrix () - Calculate Gram matrix
      %
      % Compute the Gram matrix of D
      % X = [x_1 - x_N; ... x_N-1 - x_N]
      % G = X X^t
      % G[i,j] = <x_i - x_N, x_j - x_N>
      %        = d_iN d_jN cos \beta_ij [opposite angle to d_ij]
      %        = 1/2 (d_iN^2 + d_jN^2 - dij^2)
      %
      % Schoenberg 1935: D is realizable in d-space iff Gram matrix G is positive semidefinite of rank d
      %
      % Cholesky decomposition: G is positive semidefinite of rank d iff it is Cholesky-decomposable,
      % i.e. G=XX^t where X is lower-triangular of size N-1 x d
      %
      %   obj    - Distance space object
      %   G      - Gram matrix
      %   lambda - Eigenvalues o Gram matrix (optional)
      %
      for k=1:obj.N
        for l=1:obj.N
          G(k,l) = (obj.D(k,obj.N)^2 + obj.D(l,obj.N)^2 - obj.D(k,l)^2)/2;
        end
      end
      if nargout > 1
        lambda = eig (G);
      end
    end
  end
end