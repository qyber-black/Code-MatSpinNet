function [p,q,short] = sim_dio_approx (theta,x,shortest)
  % [p,q,short] = sim_dio_approx (theta,x,shortest) - Simultaneous Diophantine Approximation
  %
  % On input, this function takes a d-dimensional COLUMN vector \theta of possibly irrational numbers
  % and on return it provides a "good," albeit not Dirichlet-good, diophantine approximation
  %
  %   theta ~ p/q
  %
  % where p is a d-dimensional vector of integers and q the integer common denominator.
  %
  % This diophantine approximation utilizes x as scaling parameters.
  % If x=s is a single number, then uniform scaling will be used. If x is a
  % (d+1)-dimensional vector [y s] then y is the non-uniform scaling vector and
  % s is a scalar scaling parameter. The dimension d=ceil((N+1)/2)-2 must be
  % entered manually as the number of components of x in this case.
  %
  % The parameter s should be decreasing to zero for better quality approximation.
  % Also on return, short is a "short vector" in the lattice
  %
  % B_s(\theta)\mathbb{Z}^{d+1}
  %
  % This vector somehow provides an idea of the approximation error.
  %
  % If shortest is given and not 0, then find shortest vector from LLL approximation
  % (for uniform scaling).
  %
  % Contrary to the paper by Kovacs that uses theta-p/q as error, here we use the error
  % theta*q-p because that is the one relevant to attainability. The returned
  % error is
  %
  % error=theta*q-p
  %
  % This fourth version of the simultaneous diophantine approximation is inspired from the paper
  %
  % A. Kovacs and N. Tihanyi,
  % "Efficient computing of n-dimensional simultaneous diophantine
  % approximation problems."
  % Acta Univ. Sapientiae, Informatica,
  % Volume 5, Number 1,
  % pp.16-34, 2013.
  %
  %*************************************************************************
  % WARNING: This paper is sloppily written.
  % While it is fundamentally correct, the way (p,q) is derived from the
  % short vector of the reduced basis of the lattice is wrong.
  % The authors forgot several scalings.
  % I had to rederive (p,q) my own way.
  %*************************************************************************
  %
  % The above paper uses a scalar parameter x that should go to infinity for
  % accrued accuracy. The significant departure from the paper is that here
  % we make the scaling a diagonal matrix X, the diagonal matrix expansion of
  % the d-dimensional vector x on input.
  % In addition, in the above paper, B(d+1,d+1)=1. Here we add the extra
  % parameter s and set B(d+1,d+1)=s. Thus if we set x=ones(d), we recover
  % the original version (attainability.m).
  %
  % Also on return, "short" is a short vector in the lattice
  %
  % B_s(\theta)\mathbb{Z}^{d+1}
  %
  % This vector somehow provides an idea of the approximation error.
  %
  if nargin < 3
    shortest = 0;
  end
  d=length(theta);
  %

  if length(x) > 1
    % Non-uniform scaling
    if length(x) ~= d + 1
      error (['x parameter must have ' d+1 ' entries or have length 1 for uniform scaling']);
    end
    X=diag(x(1:d));
    %
    Bstheta=[X -X*theta; zeros(1,d) x(d+1)];
  else
    % uniform scaling
    Bstheta=[eye(d,d) -theta; zeros(1,d) x];
  end
  %
  % Now we compute a reduced basis of the lattice
  %
  B = LLL_reduction(Bstheta);

  % Short vector
  short=B(1:end,1);

  % Find shortest vector
  if shortest ~= 0 && length(x) == 1
    D = d+1;
    bound = floor((2/sqrt(3))^D);
    len = norm(short,2);
    l=-bound*ones(D,1);
    while 1
      if sum(l) ~= 0
        X = B*l;
        Xlen = norm(X,2);
        if Xlen < len
          len = Xlen;
          short = X;
        end
      end
      n=1;
      l(1) = l(1) + 1;
      while l(n) > bound
        l(n) = -bound;
        n = n + 1;
        if n > D
          break
        end
        l(n) = l(n) + 1;
      end
      if n > D
        break
      end
    end
  end

  %
  if length(x) > 1
    % Non-uniform scaling
    q = nearest (short(d+1,1) / x(d+1));
    p = nearest ( (q * theta.*x(1:d)' + short(1:d))./x(1:d)');
  else
    % Uniform scaling
    q=nearest(short(d+1,1)/x);
    p=nearest(short(1:d,1)+q*theta);
  end
  %
  short=short(1:d,1);

end
