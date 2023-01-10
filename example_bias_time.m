function [best,results] = example_bias_time()
  % QSN package example
  %
  % Check best static control for each time and transition

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  % Construct a ring network
  N = 9;
  fprintf('=======================================================================\nConstruct an XX ring network with %d nodes\n', N);
  ring = qsn.QSN ('ring', N)

  B = NaN;
  symm = 0;

  % Setup bias control problem
  if symm ~= 0
    if init < target
      K = floor((target - init) / 2);
      for l = 0:K
        M = init+l;
        C{M} = zeros(N,N);
        C{M}(M,M) = 1;
        C{M}(target-l,target-l) = 1;
        M
        C{M}
      end
      K = ceil(N/2) - K - 1;
      for l = 1:K
        M = M + 1;
        p = init-l;
        if p < 1
          p = N+p;
        end
        C{M} = zeros(N,N);
        C{M}(p,p) = 1;
        p = target+l;
        if p > N
          p = p - N + 1;
        end
        C{M}(p,p) = 1;
        M
        C{M}
      end
      C
    else
      error('Not implemented');
    end
  else
    C = cell(1,N);
    for l=1:N
      C{l} = zeros (N,N);
      C{l}(l,l) = 1;
    end
  end

  Repeats = 10;
  Ts = 1:.2:5;
  Targets = 2:ceil(N/2);

  XX = ones(size(Targets,2),1)*Ts;
  YY = Targets'*ones(1,size(Ts,2));
  results = zeros(ceil(N/2)-1,size(Ts,2),Repeats);
  best = ones(ceil(N/2)-1,size(Ts,2));

  figure(1);
  clf;
  cm = parula(128);

  init = 1;
  for target = 1:size(Targets,2)

    subplot(ceil(size(Targets,2)/2),2,target);
    xlabel('Time');
    ylabel('Error');
    title(sprintf('Best error for 1-%d in %d-ring', Targets(target), N));
    hold on;

    for T = 1:size(Ts,2)

      ctrl = qsn.Control (ring, 'state_transfer', C,  init, Targets(target), Ts(T), B)
      for l=1:Repeats
        [x,err,exec_time,exit_flag,output,x0] = ctrl.solve_static (qsn.Opt ('fmin'), []); % 2 or 3 for Opt for fminunc, fmincon only or 8 for DE
        results(target,T,l) = err;
        if err < best(target,T)
          best(target,T) = err;
        end
        plot(Ts(T),err,'.','Color',cm(1+round(err*(size(cm,1)-1)),:));
        drawnow();
        refresh();
      end

    end

  end

end