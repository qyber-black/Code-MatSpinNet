function DD = fix_triangle_violations_mincon (D, mode, prob, alg, show)
  % DD = fix_triangle_violations_mincon (D, mode, sqp, show)
  %
  % Fix triangle violations in distance matrix D using linear programming.
  %
  %   D     - Distance matrix
  %   mode  - Determines norm for optimisation:
  %             mode > 0 - L_mode norm
  %             mode < 0 - L_mode norm with multiplictive distortion
  %             mode = [N,1] - L_N norm with scaling factor
  %   alg    - Determines algorithm for constraint optmisation:
  %              alg = 0 trust-region reflective
  %              alg = 1 interior point
  %              alg = 2 sqp
  %   prob  - If true, minimise probability change, not distance change
  %   show  - If true, verbose mode
  %   DD    - Resulting distance matrix
  %

  n = size(D,1);

  S = (n * n - n) / 2;

  Xorig = [];
  for c = 1:n-1;
    Xorig = [Xorig; D(c+1:n,c)];
  end

  X0 = ones (S, 1);

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

  B = zeros(c,1);

  LB = zeros(S, 1);
  UB = Inf(S, 1);

  if alg == 0
    alg='trust-region-reflective';
  elseif alg == 1
    alg = 'interior-point';
  else
    alg = 'sqp';
  end
  if show == 0
    dp = 'off';
  else
    dp = 'Iter';
  end
  options = optimset('Algorithm', alg, 'Display', dp);

  if prob == 0
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
  else
    if size(mode,2) == 2
      nn = length(X0);
      if mode(1) >= 0
        mode = mode(1);
        func = @(x)prob_error_scaling (x);
      else
        mode = -mode(1);
        func = @(x)prob_error_scaling_mult (x);
      end
    elseif mode >= 0
      func = @(x)prob_error (x);
    else
      mode = -mode;
      func = @(x)prob_error_mult (x);
    end
    XorigExp = exp(-Xorig);
  end

  [X,FVAL,EXITFLAG,OptDetails] = fmincon(func, X0, A, B, [], [], LB, UB, [], options);

  if show ~= 0
    display(OptDetails);
  end

  DD = zeros (n, n);
  p = 1;
  l = n-1;
  for c = 1:n-1
    DD(c+1:n,c) = X(p:p+l-1,1);
    p = p + l;
    l = l - 1;
  end
  DD = DD + DD';

  if show ~= 0
    display(sprintf('Maximum distance change: %f', max(max(abs(D-DD)))));
    display(sprintf('Norm of change to transition probability matrix: %f',norm(exp(-D)-exp(-DD))));
  end

  function e = distance_error (x)
    e = double(norm (Xorig - x, mode));
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

  function e = prob_error (x)
    e = norm (XorigExp - exp(-x), mode);
  end

  function e = prob_error_mult (x)
    e = norm ((XorigExp - exp(-x))./XorigExp, mode);
  end

  function e = prob_error_scaling (x)
    e = norm (XorigExp * x(nn) - exp(-x(1:nn-1)), mode);
  end

  function e = prob_error_scaling_mult (x)
    e = norm ((XorigExp * x(nn) - exp(-x(1:nn-1)))./XorigExp, mode);
  end

  function p = get_ind(r, c, n)
    % r > c !
    p = sum(n-c+1:n-1) + (r - c);
  end

end