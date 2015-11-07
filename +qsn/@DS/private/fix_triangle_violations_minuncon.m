function DD = fix_triangle_violations_minuncon (D, E, mode, show)
  % DD = fix_triangle_violations_minuncon (D, E, mode, show)
  %
  % Fix triangle violations in distance matrix D using optimisation.
  %
  %   D     - Distance matrix
  %   E     - Precision
  %   mode  - Determines norm for optimisatioN:
  %             mode > 0 - L_mode norm
  %             mode < 0 - L_mode norm with multiplictive distortion
  %             mode = [N,1] - L_N norm with scaling factor
  %   show  - If true, verbose mode
  %   DD    - Resulting distance matrix
  %

  n = size(D,1);

  S = (n * n - n) / 2;
  Xorig = [];
  for c = 1:n-1;
    Xorig = [Xorig; D(c+1:n,c)];
  end
  %X0 = ones (S, 1);
  X0 = double(Xorig);
  if size(mode,2) == 2
    X0 = [X0; 1];
    S = S + 1;
  end

  C = nchoosek (1:n, 3);
  cc = size(C,1);
  A = zeros (cc * 3,S);
  c = 0;
  for t = 1:cc
    %  D(b,a) + D(c,b) - D(c,a) > 0
    %  D(b,a) - D(c,b) + D(c,a) > 0
    % -D(b,a) + D(c,b) + D(c,a) > 0
    A(c + [1:3],[get_ind(C(t,2),C(t,1),n) get_ind(C(t,3),C(t,2),n) get_ind(C(t,3),C(t,1),n)]) = [ -1 -1 1; -1 1 -1; 1 -1 -1];
    c = c + 3;
  end

  if show == 0
    dp = 'off';
  else
    dp = 'Iter';
  end

  if size(mode,2) == 2
    nn = length(X0);
    if mode(1) >= 0
      mode = mode(1);
      func = @(x)distance_error_scaling (x);
    else
      mode = -mode(1);
      func = @(x)distance_error_scaling_mult (x);
    end
  elseif mode >= 0
    func = @(x)distance_error (x);
  else
    mode = -mode;
    func = @(x)distance_error_mult (x);
  end

  options = optimset('Display', dp);
  [X,FVAL,EXITFLAG,OptDetails] = fminunc(@(x)err_func(x), X0, options);

  if show ~= 0
    display(OptDetails);
  end

  DD = zeros (n, n);
  p = 1;
  l = n-1;
  for c = 1:n-1;
    DD(c+1:n,c) = X(p:p+l-1,1);
    p = p + l;
    l = l - 1;
  end
  DD = DD + DD';

  if show~= 0
    display(sprintf('Maximum distance change: %f', max(max(abs(D-DD)))));
  end

  function e = err_func (x)
    tv = A * x;
    tv(tv>=E) = 0;
    e = double(func(x) + norm(tv, mode));
  end

  function e = distance_error (x)
    e = norm (Xorig - x, mode);
  end

  function e = distance_error_mult (x)
    e = norm ((Xorig - x)./Xorig, mode);
  end

  function e = distance_error_scaling (x)
    e = norm (Xorig * x(nn) - x(1:nn-1), mode);
  end

  function e = distance_error_scaling_mult (x)
    e = norm ((Xorig * x(nn) - x(1:nn-1))./Xorig, mode);
  end

  function p = get_ind(r, c, n)
    % r > c !
    p = sum(n-c+1:n-1) + (r - c);
  end

end