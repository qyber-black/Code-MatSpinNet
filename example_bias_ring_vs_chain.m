function example_time_ring_chain
  % QSN package example
  %
  % Probabilities for information transfer from spin 1 to spin k in a
  % chain of length k, for k=1:ceil(N/2) vs 1-k transition in B-quenched
  % N-ring up to maximum time.

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  N = 13;
  Tmax = 30;
  B = 10;

  figure(1);
  clf;

  T = 0:.01:Tmax;
  Nd = ceil(N/2);
  for c = 2:Nd
    % Chain 1-c
    chain = qsn.QSN('chain', c, 'XX');
    chain_P = chain.prob(T);
    chain_1c = arrayfun(@(n) chain_P{n}(1,c),1:length(chain_P));
    % Quenched N-ring for 1-c transition
    Z = [zeros(c,1); B*ones(N-c,1)];
    ring = qsn.QSN('ring', N, 'XX', Z);
    ring_P = ring.prob(T);
    ring_1c = arrayfun(@(n) ring_P{n}(1,c),1:length(ring_P));
    % Plot
    subplot(ceil((Nd-1)/2),2,c-1);
    plot(T,ring_1c,'b');
    hold on;
    plot(T,chain_1c,'r');
    hold off;
    title(sprintf('1-%d transition in %d-chain (red) vs %g-quenched %d-ring (blue)', c, c, B, N));
  end

end