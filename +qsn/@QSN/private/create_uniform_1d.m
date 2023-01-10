function create_uniform_1d (obj, ring, in)
  % create_uniform_1d (obj, ring, in) - Create a 1D spin network
  %
  % Construct a physically uniform quantum spin network ring or chain
  % considering the first excitation sub-space only.
  %
  %   ring   - indicates ring (~= 0) or chain (== 0)
  %   in     - Cell array describing 1D network as follows:
  %            in{1} - Number of spins or vector with coupling strengths between
  %                      neighbouring spins
  %            in{2} - 'XX': XX Coupling
  %                    'H': Heisenberg coupling (optional, default 'XX')
  %            in{3} - Vector giving magnetic field strengths on individual spins
  %                    for engineered 1D networks
  %   obj    - Quantum spin network object
  %

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  % Check arguments
  if size(in,3) > 3
    error ('Too many parameters for uniform 1D spin network');
  end

  % Construct initial Hamiltonian
  if size(in,2) < 1
    error ('Number of spins or coupling strengths not defined');
  end
  d = in{1};
  if isscalar(d) && d == double(uint64(d))
    if d < 2
      error ('There must be at least two spins');
    end
    N = d;
    if ring
      d = ones (1,N);
    else
      d = ones (1,N-1);
    end
  elseif isvector(d) && ~isstr(d)
    if ring
      N = length(d);
    else
      N = length(d) + 1;
    end
    if N < 2
      error ('There must be at least two spins');
    end
  else
    error ('Illegal number of spin specification');
  end

  if ring
    H = diag(d(1:N-1),1) + diag(d(1:N-1),-1);
    H(1,N) = d(N);
    H(N,1) = d(N);
  else
    H = diag(d,1) + diag(d,-1);
  end

  % Modify Hamiltonian according to coupling type
  if size(in,2) > 1
    if ~isstr(in{2})
      error('Quantum spin network coupling type must be a string');
    end
    if strcmpi(in{2},'XX')
      ;
    elseif strcmpi(in{2},'H')
      % Heisenberg coupling in first excitation sub-space
      H = H + diag (1/2 * sum(sum(triu(H))) * ones(N,1) - sum(H,2));
    else
      error (['Unknown spin coupling type ' in{2}]);
    end
  end

  % Add engineered network fields
  if size(in,2) > 2
    Z = in{3};
    if ~isvector (Z)
      error ('Engineered network argument must be a vector');
    end
    if size(Z,1) > 1
      Z = Z';
    end
    if size(Z,2) < N
      Z = [Z ones(1,N-size(Z,2))];
    elseif size(Z,2) ~= N
      error (sprintf('Engineered ring vector length (%d) must be equal to the number of spins (%d)', size(Z,2), N));
    end
    H = H + diag(Z);
  end

  % Construct physical positions of equally spaced spin ring
  if ring
    P = 2.0 * pi * [0:N-1] / N;
    P = [cos(P); sin(P)];
  else
    P = [ ([0:N-1] / (N-1) * 2) - 1; zeros(1,N) ];
  end

  % Store data
  obj.N = N;
  obj.H = H;
  obj.pos = P;
  obj.excite = 1;

end