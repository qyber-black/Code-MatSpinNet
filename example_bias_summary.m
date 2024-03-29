function [shortest_err,shortest_time] = example_bias_summary(shortest_time)
  % function [shortest_err,shortest_time] = example_bias_summary(shortest_time)

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  min_err = 10^-2;

  for N = 5:15
    for target = 2:ceil(N/2)
      fprintf('Err: %g, N = %d, target = %d\n', min_err, N, target);
      if isnan (shortest_time(N,target))
        [shortest_err(N,target),shortest_time(N,target),~] = find_ring_bias(N,target,min_err);
      end
      fprintf('Time: %g\n', shortest_time(N,target));
    end
  end

end