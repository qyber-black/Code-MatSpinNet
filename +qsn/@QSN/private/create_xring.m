function create_xring (obj, in)
  % create_xring (obj, in) - Create an xring spin network
  %
  % Construct a physically uniform quantum spin network ring with chain
  % extensions considering the first excitation sub-space only.
  %
  %   in     - Cell array describing xring network as follows:
  %            in{1} - Number if spins in ring
  %            in{2} - Step between chains on the ring
  %            in{3} - Length of chains
  %            in{4} - 'XX': XX Coupling
  %                    'H': Heisenberg coupling (optional, default 'XX')
  %            in{5} - Vector giving magnetic field strengths on individual spins
  %                    for engineered xring networks
  %   obj    - Quantum spin network object
  %

  % Check arguments
  if size(in,2) < 3 || size(in,2) > 5
    error ('Wrong number of parameters for star spin network');
  end

  % Construct initial Hamiltonian

  d_ring = in{1};
  if isscalar(d_ring) && d_ring == double(uint64(d_ring))
    if d_ring < 2
      error ('There must be at least two spins');
    end
    N = d_ring;
    d_ring = ones (1,N);
  elseif isvector(d_ring) && ~isstr(d_ring)
    N = length(d_ring);
    if N < 2
      error ('There must be at least two spins');
    end
  else
    error ('Illegal number of ring spin specification');
  end

  M = in{2};
  if isscalar(M) && M == double(uint64(M))
    if M < 1
      error ('Cannot have less than 1 step between chains');
    end
  else
    error ('Illegal number of steps between chains');
  end

  d_chain = in{3};
  if isscalar(d_chain) && d_chain == double(uint64(d_chain))
    if d_chain < 1
      error ('There must be at least one spin in chain');
    end
    L = d_chain;
    d_chain = ones (1,L-1);
  elseif isvector(d_chain) && ~isstr(d_chain)
    L = length(d_chain) + 1;
    if L < 1
      error ('There must be at least one spin in chain');
    end
  else
    error ('Illegal number of spin in chain specification');
  end

  % Construct Hamlitonian
  S = N + ceil(N/M) * L;
  H = zeros(S,S);

  % Ring
  H(1:N,1:N) = diag(d_ring(1:N-1),1) + diag(d_ring(1:N-1),-1);
  H(1,N) = d_ring(N);
  H(N,1) = d_ring(N);

  % Chains
  p = N + 1;
  for chain = 1:M:N
    q = p + L - 1;
    H(p:q,p:q) = diag(d_chain,1) + diag(d_chain,-1);
    H(p,chain) = 1;
    H(chain,p) = 1;
    p = p + L;
  end

  % Modify Hamiltonian according to coupling type
  if size(in,2) > 3
    if ~isstr(in{4})
      error('Quantum spin network coupling type must be a string');
    end
    if strcmpi(in{4},'XX')
      ;
    elseif strcmpi(in{4},'H')
      % Heisenberg coupling in first excitation sub-space
      H = H + diag (1/2 * sum(sum(triu(H))) * ones(S,1) - sum(H,2));
    else
      error (['Unknown spin coupling type ' in{2}]);
    end
  end

  % Add engineered network fields
  if size(in,2) > 4
    Z = in{5};
    if ~isvector (Z)
      error ('Engineered network argument must be a vector');
    end
    if size(Z,1) > 1
      Z = Z';
    end
    if size(Z,2) < S
      Z = [Z ones(1,S-size(Z,2))];
    elseif size(Z,2) ~= S
      error (sprintf('Engineered ring vector length (%d) must be equal to the number of spins (%d)', size(Z,2), S));
    end
    H = H + diag(Z);
  end

  % Construct physical positions of equally spaced spin ring
  P = 2.0 * pi * [0:N-1] / N;
  P = [cos(P); sin(P)];
  for chain = 1:M:N
    P = [P, P(1:2,chain) * [2:2+L-1]];
  end

  % Store data
  obj.N = S;
  obj.H = H;
  obj.pos = P;
  obj.excite = 1;

end