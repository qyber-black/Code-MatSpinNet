classdef (Sealed = true) QSNRange < handle
  % QSNRange   Range of quantum spin networks class
  %
  % This class represents a range of QSN objects.
  %
  % QSNRange Properties:
  %   ctor      - Generator function of QSN objects given arguments from R
  %               (of the form @(x,...) QSN(...,x,...))
  %   R         - Cell array of ranges
  %   T         - Cell array of range names
  %   QSN       - Cell array of quantum spin networks according to R
  %   N         - Cell array of number of elements in each range
  %   dim       - Number of ranges
  %
  % QSN Methods:
  %   obj = QSNRange (gen, range1, name1, ...) - C'tor
  %   obj.disp ()                              - Display quantum spin network range
  %   x = obj . exec (func)                    - Execute func on all objects
  %
  properties (SetAccess = private)
    ctor = 'undef';  % Generator of QSN object
    R    = NaN;      % Ranges
    T    = NaN;      % Names
    QSN  = NaN;      % QSN objects
    N    = 0;        % Size of ranges
    dim  = 0;        % Number of ranges
  end
  properties (GetAccess = { ?qsn.DSRange }, SetAccess = private)
    L    = NaN;      % Index grid
  end
  methods
    function obj = QSNRange (gen, varargin)
      % obj . QSNRange (gen, range1, name1, range2, name2, ...) - Create quantum spin network range
      %
      %   gen      - generator function
      %   nameN    - Name of N-th range
      %   rangeN   - N-th range as vector
      %   obj      - Quantum spin network range object
      %
      obj.ctor = gen;
      obj.T = { varargin{1:2:end} };
      obj.R = { varargin{2:2:end} };
      obj.N = cellfun(@length, obj.R);
      obj.dim = length(obj.N);
      if length(obj.N) == 1
        obj.N = [obj.N 1];
      end
    end
    function disp (obj)
      % obj.disp () - Display quantum spin network range
      %
      %   obj - Quantum spin network range object
      %
      fprintf (1, 'Quantum Spin Network Range over %d ranges\n\n Constructor: ', obj.dim);
      disp(obj.ctor);
      str = sprintf ('  %%%ds - %%d elements from %%d to %%d\n', max(cellfun(@(x) size(x,2), obj.T)) + 1);
      for l = 1:obj.dim
        fprintf (1, str, obj.T{l}, obj.N(l), obj.R{l}(1), obj.R{l}(end));
      end
      fprintf (1, '\n Grid');
      if (~iscell(obj.L))
        fprintf (1, ': not generated\n\n');
      else
        str = strtrim (sprintf('%d ', size(obj.L)));
        fprintf (1, ' size:    [%s]\n', str);
      end
      fprintf (1, ' Object');
      if (~iscell(obj.QSN))
         fprintf (1, ': not generated\n\n');
      else
        str = strtrim (sprintf('%d ', size(obj.QSN)));
        fprintf (1, ' size: [%s]\n', str);
      end
      fprintf (1, '\n');
    end
    function x = exec (obj, func)
      % x = obj . exec (func) - Execute func on all objects
      %
      % Executes func of the form @x x.function (args...) on all objects in the range
      % and returns the result in a cell array x (of the same form than the R cell array).
      %
      %   func - Function handle of form @x x.function (args...) to execute on objects x
      %   obj  - Object of DSRange class
      %   x    - Cell array of function retuns (optional)
      %
      if ~isa(func, 'function_handle')
        error ('Illegal function handle');
      end
      if nargout
        x = cellfun(func, obj.QSN);
      else
        cellfun(func, obj.QSN);
      end
    end
  end
  methods (Access = { ?qsn.DSRange })
    function generate_grid (obj)
      % obj . generate_grid () - Generate a grid from the given ranges (and cache results)
      %
      %   obj - Object of QSNRange class
      if ~iscell(obj.L)
        obj.L = cell(obj.N);
        l = ones(obj.dim, 1); % Initialise dim-D counters for ranges (dim-D index)
        ll = 0;               % Current counter (sequential index of dim-D index)
        % Simulate dim-D nested loops over dim-D index by increasing sequential index
        while l(obj.dim) <= obj.N(obj.dim)
          ll = ll + 1;
          obj.L(ll) = { l' }; % ll-th (sequential index) element of grid object is dim-D index l
          % Increase dim-D index l according to maximal indices in ranges
          k = 1;
          l(k) = l(k) + 1;
          while k < obj.dim && l(k) > obj.N(k)
            l(k) = 1;
            k = k + 1;
            l(k) = l(k) + 1;
          end
        end
      end
    end
    function generate_qsn(obj)
      % obj . generate_qsn ()  - Generate QSN objectes from dim-D range index (and cache results)
      %
      %   obj - Object of QSNRange class
      %
      if ~iscell(obj.QSN)
        obj.QSN = cell(obj.N);
        cellfun(@(x) obj.generate_qsn_obj(x),obj.L);
      end
    end
  end
  methods (Access = private)
    function generate_qsn_obj (obj, l)
      % obj . generate_qsn_obj ()  - Convert dim-D range index into generator argument and create QSN object
      %
      %   obj - Object of QSNRange class
      %
      P = num2cell(l);
      X = num2cell(arrayfun(@(x,y) obj.R{x}(y), 1:obj.dim, l));
      obj.QSN(sub2ind(obj.N,P{:})) = { obj.ctor(X{:}) };
    end
  end
end