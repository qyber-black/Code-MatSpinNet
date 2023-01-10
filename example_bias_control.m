function [Best_X,Best_T,Fastest_X,Fastest_T] = example_bias_control()
  % QSN package example
  %
  % Static bias control for information transfer from spin 1 to spin target
  % in N-ring either for fixed time T (T ~= NaN) and with maximum bias
  % constraint (B ~= NaN). Controls can be arbitrary or symmetric w.r.t.
  % smmetry axis between 1-target.

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  N = 11;              % Ring size
  xring_length = 0;    % Extend each spin in ring by xring_length chain (or 0 for ring only)
  ring_target = 5;     % 1 to target transition (in ring!)

  T = NaN;             % Fix time?
  B = NaN;             % Maximum bias?
  readout = [0 0 400]; % Readout window (or empty to ignore)
                       % readout(1) Time window for readout (0 - do not consider window)
                       % readout(2) Pertubation (or 0 for none, time window required)
                       % readout(3) Samples for expectation over pertubation (if pertubation not 0)
  min_err = 0.01;      % Largest error acceptable for shortest solution.
  repeats = 100;       % How many restarts?

  symm = 1;            % Enforce symmetry?
  initT = [30 500];    % (1) - max time, (2) - samples for init times from chain (or empty)
  bias_init = 1;       %  Initial biases with peak/trough/constant instead of random (0)?
  noise = 0;           % Add noise to initial values(or 0)

  % Setup network
  if xring_length == 0
    ring = qsn.QSN ('ring', N);
    target = ring_target;
  else
    ring = qsn.QSN ('xring', N, 1, xring_length);
    target = N + xring_length * ring_target;
  end

  % Find bias controls
  [Results,best,fastest] = ring.bias_control(1,target, T, B, readout, min_err, repeats, symm, initT, bias_init, noise);

end