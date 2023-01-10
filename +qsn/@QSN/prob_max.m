function [pm,tm] = prob_max (obj, tmin, tmax, steps, depth, show)
  % [pm,tm] = obj . prob_max (tmin, tmax, steps, depth, show)
  %  - Estimate maximum transition probability in time interval
  %
  %   tmin   - Start time
  %   tmax   - End time
  %   steps  - Number of steps to subdivide time interval
  %   depth  - Recursion depth
  %   show   - 0: do not display,
  %            1: display final result,
  %            2: display updates and final result
  %   obj    - Quantum spin network object
  %   pm     - Estimate of maximum transition probability
  %   tm     - Time when estimate is attained

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  if ~isscalar(tmin) || tmin < 0 || isnan(tmin) || isinf(tmin)
    error ('Illegal tmin value');
  end
  if ~isscalar(tmax) || tmax < tmax || isnan(tmax) || isinf(tmax)
    error ('Illegal tmax value');
  end
  if ~isscalar(steps) || steps < 1 || steps ~= double(uint64(steps))
    error ('Illegal steps value');
  end
  if ~isscalar(depth) || depth < 1 || depth ~= double(uint64(depth))
    error ('Illegal depth value');
  end
  if ~isscalar(show) || (show ~= 0 && show ~= 1 && show ~= 2)
    error ('Ilelgal show value');
  end
  
  % Initialise output
  pm = -Inf(obj.N);
  tm = Inf(obj.N);

  % Initialise helpers
  i = zeros(obj.N);

  % Check interval
  check (tmin, tmax, depth);

  % Render final result
  if show ~= 0
    render (Inf);
  end

  function check (t0, t1, d)
    % Recursively check interval [t0,t1] for maximum at recursion level d.

    % Step size and terminate if step size is too small
    h = (t1-t0) / steps;
    h2 = h / 2;
    if h < obj.eps
      warning (sprintf ('Terminating recursive search early:\n Step size %f\n Time interval [%f,%d]\n Recursion depth: depth-d', h, t0, t1));
      return
    end

    % Check t0
    update (t0, 0, h2, d);

    % Check probabilities between t0 and t1
    for t = [t0+h:h:t1-h]
      update (t, h2, h2, d);
    end

    % Check t1
    update (t1, h2, 0, d);
  end

  function r = update (t, d0, d1, d)
    % Update probabilities from values at t and check interval [t-d0,t+d1] if change occured
    p = obj.prob (t);
    % Check if update of maxima is necessary
    i = p>pm;
    if max(max(i))
      pm(i) = p(i);
      tm(i) = t;
      if show > 1
        render(t)
      end
      % Explore recursively as update occured
      if d >= 0
        check (t-d0, t+d1, d-1);
      end
    end
  end

  function render(t)
    % Render the maxima with their times
    figure (1);
    clf;
    subplot(2,1,1);
    imagesc(pm);
    if isinf(t)
      title(sprintf('Maximum transition probabilities in [%15.5f,%15.5f]', tmin, tmax));
    else
      title(sprintf('Maximum transition probabilities in [%15.5f,%15.5f] (updated at %15.5f)', tmin, tmax, t));
    end
    colorbar;
    subplot(2,1,2);
    imagesc(tm);
    title('Time when maximum transition probability is attained');
    colorbar;
    pause (.01);
  end

end