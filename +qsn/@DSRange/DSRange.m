classdef (Sealed = true) DSRange < handle
  % DSRange   Range of distance spaces class
  %
  % This class represents a range of DS objects
  %
  % DSRange Properties:
  %   R         - Cell array of ranges
  %   T         - Cell array of range names
  %   DS        - Cell array of distance spaces according to R
  %   N         - Cell array of number of elements in each range
  %   dim       - Number of ranges
  %
  % DSRange Methods:
  %   obj = DSRange (QSNR)                         - Create a distance space range from a DSRange
  %   [delta4,scaled_delta4] = obj . max_Gromov4pt (verbose, d, l)
  %                                                - Maximum Gromov deltas over range
  %   [max_violation] = obj . max_triangle_inequality_violation (verbose, d, l)
  %   obj . plot_triangle_inequality_violation (l) - Plot triangle inequality violations
  %   x = obj . exec (func)                        - Execute func on all objects
  %

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  properties (SetAccess = private)
    R    = NaN;      % Ranges
    T    = NaN;      % Names
    N    = 0;        % Size of ranges
    dim  = 0;        % Number of ranges
    DS   = NaN;      % DS objects
  end
  properties (Access = private)
    L    = NaN;      % Index grid
  end
  methods
    function obj = DSRange (QSNR)
      % obj = DSRange (QSNR) - Create a distance space range from a DSRange
      %
      % This c'tor is intended to be used only by the QSNRange object to generate all
      % the distances spaces associated with the QSNs.
      %
      %   QSNR     - QSNRange object
      %   obj      = DSRange object
      %
      if ~isa(QSNR, 'qsn.QSNRange')
        error ('First argument must be a qsn.QSNRange object');
      end
      obj.R = QSNR.R;
      obj.T = QSNR.T;
      obj.N = cellfun(@length, obj.R);
      obj.dim = length(obj.N);
      if length(obj.N) == 1
        obj.N = [obj.N 1];
      end
      QSNR.generate_grid ();
      obj.L = QSNR.L;
      QSNR.generate_qsn ();
      obj.generate_ds (QSNR.QSN);
    end
    function [delta4,scaled_delta4] = max_Gromov4pt (obj, verbose, d, l)
      % [delta4,scaled_delta4] = obj . max_Gromov4pt (verbose, d, l) - Maximum Gromov deltas over range
      %
      % Calculated the 4-point Gromov delta and the scaled 4-point Gromov delta over the
      % QSNs in the range and also enables plotting for 1D and 2D ranges.
      %
      %   verbose       - If true, plot Gromov deltas (optional, default: false)
      %   d             - Multiple curves (true) or 2D plot (false) for 2D ranges (optional, default: false)
      %   l             - Print legend for 2D curve plots (optional, default: false)
      %   obj           - DSRange object
      %   delta4        - Maximum Gromov 4-point delta (as in R)
      %   scaled_delta4 - Maximum scaled Gromov 4-point delta (as in R)
      %
      [delta4,scaled_delta4] = cellfun(@(x) x.max_Gromov4pt(), obj.DS, 'UniformOutput', false);
      delta4 = cell2mat(delta4);
      scaled_delta4 = cell2mat(scaled_delta4);
      if exist('verbose', 'var') && verbose
        clf;
        if obj.dim == 1
          colormap('default');
          subplot(1,2,1);
          plot(obj.R{1},delta4);
          axis([obj.R{1}(1) obj.R{1}(end) 0 max(max(delta4))]);
          xlabel(obj.T{1});
          ylabel('Maximum $\delta_4$', 'interpreter', 'latex');
          subplot(1,2,2);
          plot(obj.R{1},scaled_delta4);
          axis([obj.R{1}(1) obj.R{1}(end) 0 max(max(scaled_delta4))]);
          xlabel(obj.T{1});
          ylabel('Maximum scaled $\delta_4$', 'interpreter', 'latex');
        elseif obj.dim == 2
          if exist ('d', 'var') && d
            colormap (varycolor(obj.N(1)));
            subplot(1,2,1);
            plot (obj.R{1},delta4);
            axis([obj.R{1}(1) obj.R{1}(end) 0 max(max(delta4))]);
            xlabel(obj.T{1});
            ylabel('Maximum $\delta_4$', 'interpreter', 'latex');
            if exist ('l', 'var') && l
              L = arrayfun (@(x) sprintf('%s=%g', obj.T{2}, x), obj.R{2}', 'UniformOutput', false);
              legend(L{:});
            end
            subplot(1,2,2);
            plot(obj.R{1},scaled_delta4);
            axis([obj.R{1}(1) obj.R{1}(end) 0 max(max(scaled_delta4))]);
            xlabel(obj.T{1});
            ylabel('Maximum scaled $\delta_4$', 'interpreter', 'latex');
            if exist ('L', 'var')
              legend(L);
            end
          else
            colormap('default');
            subplot(1,2,1);
            surf(obj.R{2},obj.R{1},delta4);
            hold on;
            surf(obj.R{2},obj.R{1},zeros(size(delta4)),delta4);
            hold off;
            axis([obj.R{2}(1) obj.R{2}(end) obj.R{1}(1) obj.R{1}(end) 0 max(max(delta4))]);
            xlabel(obj.T{2});
            ylabel(obj.T{1});
            zlabel('Maximum $\delta_4$', 'interpreter', 'latex');
            subplot(1,2,2);
            surf(obj.R{2},obj.R{1},scaled_delta4);
            hold on;
            surf(obj.R{2},obj.R{1},zeros(size(scaled_delta4)),scaled_delta4);
            axis([obj.R{2}(1) obj.R{2}(end) obj.R{1}(1) obj.R{1}(end) 0 max(max(scaled_delta4))]);
            hold off;
            xlabel(obj.T{2});
            ylabel(obj.T{1});
            zlabel('Maximum scaled $\delta_4$', 'interpreter', 'latex');
          end
        else
          error('Too many ranges');
        end
      end
    end
    function [max_violation] = max_triangle_inequality_violation (obj, verbose, d, l)
      % [max_violation] = obj . max_triangle_inequality_violation (verbose, d, l)
      %
      % Calculate maximum triangle violation for all QSNs in the range.
      %
      %   verbose       - If true, plot results (optional, default: false)
      %   d             - Multiple curves (true) or 2D plot (false) for 2D ranges (optional, default: false)
      %   l             - Print legend for 2D curve plots (optional, default: false)
      %   max_violation - Maximum triangle violation (as in R)
      %
      [~,max_violation] = cellfun(@(x) x.check_triangle_inequality (), obj.DS);
      max_violation = -max_violation;
      if exist('verbose', 'var') && verbose
        clf;
        if obj.dim == 1
          colormap('default');
          plot(obj.R{1},max_violation);
          axis([obj.R{1}(1) obj.R{1}(end) 0 max(max(max_violation))]);
          xlabel(obj.T{1});
          ylabel('Maximum triangle inequality violation');
        elseif obj.dim == 2
          if exist ('d', 'var') && d
            colormap (varycolor(obj.N(1)));
            plot (obj.R{1},max_violation);
            axis([obj.R{1}(1) obj.R{1}(end) 0 max(max(max_violation))]);
            xlabel(obj.T{1});
            ylabel('Maximum triangle inequality violation');
            if exist ('l', 'var') && l
              L = arrayfun (@(x) sprintf('%s=%g', obj.T{2}, x), obj.R{2}', 'UniformOutput', false);
              legend(L{:});
            end
          else
            colormap('default');
            surf(obj.R{2},obj.R{1},max_violation);
            hold on;
            surf(obj.R{2},obj.R{1},zeros(size(max_violation)),max_violation);
            hold off;
            axis([obj.R{2}(1) obj.R{2}(end) obj.R{1}(1) obj.R{1}(end) 0 max(max(max_violation))]);
            xlabel(obj.T{2});
            ylabel(obj.T{1});
            zlabel('Maximum triangle inequality violation');
          end
        else
          error('Too many ranges');
        end
      end
    end
    function plot_triangle_inequality_violation (obj, l)
      % obj . plot_triangle_inequality_violation (l) - Plot triangle inequality violations
      %
      %   l - Print legend for 1D curve plots (optional, default: false)
      %
      clf;
      hold all;
      set (gca, 'ColorOrder', varycolor (prod(obj.N)));
      set(0,'defaultAxesFontName', 'Arial');
      set(0,'defaultAxesFontSize', 28);
      set(0,'defaultTextFontName', 'Arial');
      if obj.dim == 1
        for l = 1:length(obj.DS);
          [~,~,T] = obj.DS{l}.check_triangle_inequality ();
          err = sort(-reshape(T.*(T<-obj.DS{l}.epsilon ()), prod(size(T)), 1), 1, 'descend');
          L = min (find(err < obj.DS{l}.epsilon ()));
          plot ([1:L],sort (err(1:L), 1, 'descend'), 'LineWidth', 2);
        end
        title ('Triangle violations');
        xlabel ('Triangles sorted by violation');
        ylabel ('Triangle violation');
        if exist ('l', 'var') && l
          L = arrayfun (@(x) sprintf('%s=%g', obj.T{1}, x), obj.R{1}', 'UniformOutput', false);
          legend(L{:});
        end
        hold off;
      else
        error('Too many ranges');
      end
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
        x = cellfun(func, obj.DS);
      else
        cellfun(func, obj.DS);
      end
    end
  end
  methods (Access=private)
    function generate_ds(obj, QSNR)
      % obj.generate_ds (QSNR) - Generate distance spaces for QSNs 
      %
      %   QSNR - QSNRange object
      %   obj  - DSRange object
      %
      if ~iscell(obj.DS)
        obj.DS = cell(obj.N);
        cellfun(@(x,y) obj.generate_ds_helper(x,y), QSNR, obj.L);
      end
    end
    function generate_ds_helper (obj, n, l)
      % obj.generate_ds_helper (QSNR) - Generate distance space from QSN object
      %
      %   QSN - QSN object to genereate distance space
      %   l   - Index for generated DS object in obj.DS
      %   obj - DSRange object
      %
      P = num2cell(l);
      obj.DS{P{:}} = qsn.DS (n);
    end
  end
end