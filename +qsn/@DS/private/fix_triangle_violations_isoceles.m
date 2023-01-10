function DD = fix_triangle_violations_isoceles (D, eps, mode, show)
  % DD = fix_triangle_violations_isoceles (D, E, mode, show)
  %
  % Fix triangle violations in distance matrix D by taking
  % 'mode' (mean, max) of two largest violation distances
  % iteratively,
  %
  %   D    - Distance matrix
  %   eps  - Precision
  %   mode - Mode of correcting distances (@(x) mean(x), @(x) max(x))
  %   show - If true, verbose mode
  %   DD   - Resulting distance matrix
  %

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  if show ~= 0
    D0 = D;
  end
  DD = D;
  N = size(D,1);
  T = zeros(N);

  while 1,

    % Check if distances in D fulfill triangle inequality
    for k=1:N
      for l=1:N
        for j=1:N
          T(k,l,j) = D(k,l)+D(l,j)-D(j,k);
        end
      end
    end
    x = -min(min(min(T)));

    if show ~= 0
      count = 0;
      display(sprintf('Max triangle inequality violation: %d', max(x,0)));
    end

    if x < eps
      if show ~= 0
        display(sprintf('Maximum distance change: %f', max(max(abs(D0-D)))));
      end
      break
    end;

    [A B C] = ind2sub(size(T), find (T < -eps));
    X = sort([A B C],2);

    for I = 1:size(X,1)

      a = X(I,1);
      b = X(I,2);
      c = X(I,3);
      ok = 1;
      for J = 1:I-1
        if a == X(J,1) && b == X(J,2) && c == X(J,3)
          ok = 0;
          break;
        end
      end

      if ok == 1
        % [d,i] = min([D(a,b) D(b,c) D(c,a)]); % preserves smallest edge
        [d,i] = min(abs([D(b,c)-D(c,a), D(a,b)-D(c,a), D(a,b)-D(b,c)])); % preserves edge with most different length
        if i == 1 % a,b
          avg = sqrt (mode([D(b,c),D(c,a)])^2 + D(a,b)^2/2);
          DD(b,c) = avg;
          DD(c,b) = avg;
          DD(c,a) = avg;
          DD(a,c) = avg;
        elseif i == 2 % b,c
          avg = sqrt (mode([D(a,b),D(c,a)])^2 + D(b,c)^2/2);
          DD(a,b) = avg;
          DD(b,a) = avg;
          DD(c,a) = avg;
          DD(a,c) = avg;
        else % c,a
          avg = sqrt (mode([D(a,b),D(b,c)])^2 + D(c,a)^2/2);
          DD(a,b) = avg;
          DD(b,a) = avg;
          DD(b,c) = avg;
          DD(c,b) = avg;
        end
        if show ~= 0
          count = count + 1;
          e = min([T(a,b,c), T(a,c,b),T(b,a,c), T(b,c,a),T(c,a,b), T(c,b,a)]);
          display(sprintf('%3d: triangle %d %d %d / %f %f %f : %f', count, a, b, c, D(a,b), D(b,c), D(a,c), e));
          display(sprintf('     => %f %f %f', DD(a,b), DD(b,c), DD(a,c)));
        end
      end

    end

    D = DD;

  end

end
