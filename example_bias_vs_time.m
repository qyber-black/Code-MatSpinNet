function [Time,Prob,Err] = example_bias_vs_time()
  % QSN package example
  %
  % Time for max. probability given bias for ring

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  N = 13;                  % Ring size
  B = [0 10.^[1:6]];       % Bias range
  Tmax = 100000;           % Max time
  Tstep = 0.001;           % Time step
  Trange = 100;            % Length of interval for step-wise check until Tmax
  target = 6;              % Target
  bias_pos = 10;           % Position for bias

  Time = [];
  Prob = [];
  Err = [];

  n = 1;
  for k = 1:size(B,2)
    bias = zeros(1,N);
    bias(1,bias_pos) = B(k);
    ring = qsn.QSN('ring',N,'XX',bias)

    p_max = ring.prob();

    tr_prob = [];
    tr_time = [];
    t = -1;
    p = -1;
    start_time = 0;
    while start_time < Tmax
      range_time = start_time:Tstep:start_time+Trange;

      tr = ring.trace(1,target,range_time);
      [pp,ind] = max(tr);
      tr_prob = [tr_prob pp pp];
      tr_time = [tr_time start_time start_time+Trange];
      if pp > p
        p = pp;
        t = range_time(ind);
      end

      subplot(2,2,4);
      plot(tr_time,tr_prob,'r');
      title(sprintf('Current trace envelope, bias: %g',B(k)));
      pause(0.1);

      start_time = start_time + Trange;
    end

    Time(k) = t;
    Prob(k) = p;
    Err(k) = abs(p_max(1,target) - p);

    lB = log10(B(1:k));
    lB(1) = 0;

    subplot(2,2,1);
    plot(lB,Time(1:k),'g');
    title('Time vs. log(Bias)');

    subplot(2,2,2);
    plot(lB,Prob(1:k),'g');
    title('Probability vs. log(Bias)');

    subplot(2,2,3);
    plot(lB,Err(1:k),'r');
    title('|p-max - p| vs. log(Bias)');

    pause(0.1);

  end

end