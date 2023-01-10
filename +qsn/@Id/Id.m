classdef (Sealed = true) Id < handle
  % Id   Quantum spin network identification
  %
  % This class identifies quantum spin networks.
  %
  % Id Properties:
  %   S          - spin network for identification
  %   Range      - discrete parameter ranges
  %   Samples    - 1D continuous parameter space samples for discrete parameters
  %   Values     - PDF values at Samples
  %   dens       - 1D sampling density
  %   para       - Parameters for maximum
  %   paraN      - 1D para sample index for maximum
  %
  % QSN Methods:
  %   obj = Id (s, range, s_d)     - C'tor
  %   [para,err,data] = obj.find (s, Time, Exp_N, Exp_R, MaxIt)
  %                                - Find ring Hamiltonian by iterative refinement
  
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  properties (SetAccess = private)
    S          = 'undef';    % Spin network for identification
    Range      = 'undef';    % Parameter range
    Samples    = 'undef';    % Parameter space samples
    Values     = 'undef';    % PDF values at parameter space sample points
    dens       = nan;        % 1D sampling density
    para       = [nan, nan]; % Parameters for maximum
    paraN      = nan;        % 1D H sample index for maximum
  end
  properties (Access = private)
    eps        = double(eps('single')); % Precision
  end
  methods
    function obj = Id (s, range, s_d)
      % obj = Id (s) - Crate a quantum spin network identificaiton problem
      %
      %   s     - spin network for identification
      %   range - paramter range [min_N max_N; min_d max_d]
      %   s_d   - coupling strength d sample density (points per unit)
      %   obj   - Id object
      %
      if ~isa (s, 'qsn.QSN')
        error ('Illegal system to identify');
      end
      obj.S = s;
      obj.Range = cell(1);
      obj.Range{1} = range(1,1):1:range(1,2);
      obj.Samples = cell(size(obj.Range{1}));
      obj.Values = cell(size(obj.Range{1}));
      for l = 1:size(obj.Range{1},2)
        obj.Samples{l} = range(2,1):(range(2,2) - range(2,1))/(s_d-1):range(2,2);
        obj.Values{l} = zeros(size(obj.Samples{l}));
      end
      obj.dens = s_d;
    end
    function [para,err,data] = find (obj, s, Time, Exp_N, Exp_R, MaxIt)
      % [para,err,data] = obj.find (s, Time, Exp_N, Exp_R, MaxIt)
      %
      % Find ring Hamiltonian by iterative refinement.
      % Identifies length and coupling strength under assumption that
      % it is a ring with uniform coupling strength (for now - FIXME).
      %
      %    s                - spin to use for identification
      %    Time             - Time interval [T0 T1]
      %    Exp_N            - Number of experiments per time step and chain length
      %    Exp_R            - Number of times each experiment Exp_N is repeated
      %    MaxIt            - Number of iterations
      %    obj              - Id object
      %    para             - Found parameters [length coupling_strength]
      %    err              - L_0 parameter error
      %    data             - Experimental Data
      %
      err = [];   % error
      data = [];  % experimental data
      seed = [0]; % sampling seeds
      % Collect experiments and create PDF
      for it = 1:MaxIt
        disp(sprintf('Iteration %d', it));
        disp(sprintf('  New experiments: %d on trace %d-%d', Exp_N, s, s));
        % FIXME: Adative experiments? (probably not much better)
        [d,seed] = experiments(obj.S, Exp_N, Exp_R, Time, s, seed);
        data = [data, d];
        disp('  Updating PDF');
        obj.update_PDF (s, d, data);
        figure(1);
        clf;
        subplot (2,1,1);
        rn = size(obj.Range{1},2);
        sn = size(obj.Samples{1},2);
        scatter (cell2mat(obj.Samples), reshape(ones(sn,1)*obj.Range{1},1,sn*rn),20,cell2mat(obj.Values),'fill');
        colorbar;
        subplot(2,1,2);
        plot (obj.Samples{obj.paraN}, obj.Values{obj.paraN},'r*-');
        err = [err sum(abs(obj.para - [obj.S.N, obj.S.H(1,2)]))];
        disp(sprintf('    Error: %f', err(end)));
        disp('    Max. at');
        disp(obj.para);
        title (sprintf ('PDF at iteration %d with optimal N = %d',it, obj.para(1,1)));
        refresh (1);
        pause(1);
      end
      % Optimise over coupling strength
      disp('Optimise d for N');
      [obj.para(1,2),fval,ef,output] = fminsearch (@(d) -prob_binomial_log (s, [obj.para(1,1) d], data), obj.para(1,2));
      disp(output);
      err = [err sum(abs(obj.para - [obj.S.N, obj.S.H(1,2)]))];
      disp(sprintf('    Error: %f', err(end)));
      disp('    Max. at');
      disp(obj.para);
      para = obj.para;
    end
  end
  methods (Access = private)
    function update_PDF (obj, s, exp, all_exp)
      % obj.update_PDF (s, exp, all_exp) - Update probability samples
      %
      %   s       - spin used for measurement
      %   exp     - new experimental data
      %   all_exp - all experimental data
      %
      max_val = -1;
      for Ni = 1:size(obj.Range{1},2)
        N = obj.Range{1}(Ni);
        % Update probabilities with new experiments
        result = [];
        for di = 1:size(obj.Samples{Ni},2)
          d = obj.Samples{Ni}(di);
          result(1,di) = prob_binomial_log (s, [N d], exp);
        end
        obj.Values{Ni} = obj.Values{Ni} + result;
        % Resample
        ns = size(obj.Samples{Ni},2) - 2;     % number of samples inside (excl. start and end)
        V = obj.Values{Ni};                   % Values
        el = min (V);
        V = V - el;
        V = V / max(V);
        W = [0 obj.Samples{Ni}(2:end) - obj.Samples{Ni}(1:end-1)]; % Weights according to distance between samples
        X = [V(1)*W(1), V(2:end-1) .* W(2:end-1) + V(1:end-2) .* W(1:end-2), V(end)*W(end)] / sum(W);
        csum = cumsum (X);                    % cumulative sum over weighted values
        csum = csum / csum(end);              % Normalise
        csum = round (csum * ns);             % Integer steps
        pts = diff(csum);                     % number of samples per interval
        pts(end) = pts(end) + 1;              % ensure final sample is maximum of range
        new_samples = [obj.Samples{Ni}(1)];   % start at first of range
        prev = 1;
        for k = 1:size(pts,2)
          if pts(k) > 0       % place samples
            xx = obj.Samples{Ni}(k+1) - [pts(k)-1:-1:0] * (obj.Samples{Ni}(k+1) - obj.Samples{Ni}(prev)) / pts(k);
            new_samples = [new_samples xx];
            prev = k+1;
          end
        end
        obj.Samples{Ni} = new_samples;
        result = [];
        for di = 1:size(obj.Samples{Ni},2)
          d = obj.Samples{Ni}(di);
          result(1,di) = prob_binomial_log (s, [N d], all_exp);
        end
        obj.Values{Ni} = result;
        % Find maximum
        [m idx] = max(obj.Values{Ni});
        if m > max_val
          max_val = m;
          obj.paraN = Ni;
          obj.para(1,1) = N;
          obj.para(1,2) = obj.Samples{Ni}(idx);
        end
      end
    end
  end
end
