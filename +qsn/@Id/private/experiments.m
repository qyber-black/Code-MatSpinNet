function [D, seed] = experiments (S, Exp_N, Repeat_N, Time, s, seed);
  % Compute results of experiments
  %
  %    [D, seed] = experiments (S, Exp_N, Repeat_N, Time, traces, varargin);
  %
  %    S        - QSN object
  %    Exp_N    - Number of experiments per trace
  %    Repeat_N - Number of repetitions of each experiment
  %    Time     - Time interval
  %    s        - measure treace s-s
  %    seed     - Sampling seed
  %
  %    D        - Experimental data [Time...; Ones...; Zeros...]
  %    seed     - Index of sampling sequence of where to continue

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  % Sample times
  D(1,1:Exp_N) = Time(1) + i4_to_hammersley_sequence(1, Exp_N, [0], seed, [1], [2]) * (Time(2) - Time(1));
  % Execute experiments
  D(2,1:Exp_N) = sum(rand(Repeat_N, Exp_N) < ones(Repeat_N,1) * S.trace (s, s, D(1,1:end)), 1);
  D(3,1:Exp_N) = Repeat_N - D(2,1:Exp_N);
  % Increase seeds
  seed = seed + Exp_N;
end
