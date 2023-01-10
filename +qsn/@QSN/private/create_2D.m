function create_2D (obj, in)
  % create_2D (obj, in) - Create a 2D spin network
  %
  % Construct a 2D quantum spin network considering the first excitation
  % sub-space only.
  %
  %   in     - Cell array describing 2D network as follows:
  %     in{1}: N - Number if spins in ring
  %     in{2}: string == - 'random', 'xring', 'chain', 'ring' topology
  %              topology = 'random' (N is total number of spins)
  %                in{3} : D - [xmin xmax; ymin max] rectangle for positions (optional, default: [0 1; 0 1])
  %                in{4} : C - 'XX': XX Coupling
  %                            'H': Heisenberg coupling (optional, default 'XX')
  %                in{5} : L - Couplings between neighbours only (optional, default 0/false)
  %              topology = 'xring' (N is number of spins in ring)
  %                in{3} : M - Step between spins in ring that have chains
  %                in{4} : L - Length of chains
  %                in{5} : C - 'XX': XX Coupling
  %                            'H': Heisenberg coupling (optional, default 'XX')
  %                in{6} : L - Couplings between neighbours only (optional, default 0/false)
  %                in{7} : E - Vector for engineered ring (optional)
  %              topology == 'chain' (N is number of spins in chain)
  %                in{3} : C - 'XX': XX Coupling
  %                            'H': Heisenberg coupling (optional, default 'XX')
  %                in{4} : L - Couplings between neighbours only (optional, default 0/false)
  %                in{5} : E - Vector for engineered chain (optional)
  %              topology == 'ring' (N is number of spins in ring)
  %                in{3} : C - 'XX': XX Coupling
  %                            'H': Heisenberg coupling (optional, default 'XX')
  %                in{4} : L - Couplings between neighbours only (optional, default 0/false)
  %                in{5} : E - Vector for engineered ring (optional)
  %            matrix - matrix of N given positions [x0, y0; ...; xn, yn];
  %                in{3} : C - 'XX': XX Coupling
  %                            'H': Heisenberg coupling (optional, default 'XX')
  %                in{4} : L - Couplings between neighbours only (optional, default 0/false)
  %   obj    - Quantum spin network object
  %
  % FIXME: low-discrepancy positions

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  % Check arguments
  if size(in,2) < 2
    error ('Wrong number of parameters for 2D spin network');
  end

  % Number of spins
  N = in{1};
  if ~isscalar(N) || N ~= double(uint64(N)) || N < 2
    error ('Illegal number of spins');
  end

  % Construct positions
  C = 'XX';
  E = NaN;
  local = 0;
  if isstr (in{2})
    if strcmpi(in{2},'random')
      %% Generate random positions
      if size(in,2) == 5
        D = in{3};
        C = in{4};
        local = in{5};
      elseif size(in,2) == 4
        D = in{3};
        C = in{4};
      elseif size(in,2) == 3
        D = in{3};
      elseif size(in,2) == 2
        D = [0.0 1.0; 0.0 1.0];
      else
        error ('Illegal parameters');
      end
      P = [D(1,1) + (D(1,2) - D(1,1)) * rand(1,N); D(2,1) + (D(2,2) - D (2,1)) * rand(1,N)];
    elseif strcmpi(in{2},'xring')
      if size (in,2) == 4
        M = in{3};
        L = in{4};
      elseif size (in,2) == 5
        M = in{3};
        L = in{4};
        C = in{5};
      elseif size (in,2) == 6
        M = in{3};
        L = in{4};
        C = in{5};
        local = in{6};
      elseif size (in,2) == 7
        M = in{3};
        L = in{4};
        C = in{5};
        local = in{6};
        E = in{7};
      else
        error ('Illegal Paramaters');
      end
      P = 2.0 * pi * [0:N-1] / N;
      P = [cos(P); sin(P)];
      for chain = 1:M:N
        P = [P, P(1:2,chain) * [2:2+L-1]];
      end
      if size(in,2) > 4
      end
      N = N + double(uint64(N/M)) * L;
    elseif strcmpi(in{2},'chain')
      P = [ ([0:N-1] / (N-1) * 2) - 1; zeros(1,N) ];
      if size (in,2) == 3
        C = in{3};
      elseif size (in,2) == 4
        C = in{3};
        local = in{4};
      elseif size (in,2) == 5
        C = in{3};
        local = in{4};
        E = in{5};
      end
    elseif strcmpi(in{2},'ring')
      P = 2.0 * pi * [0:N-1] / N;
      P = [cos(P); sin(P)];
      if size (in,2) == 3
        C = in{3};
      elseif size (in,2) == 4
        C = in{3};
        local = in{4};
      elseif size (in,2) == 5
        C = in{3};
        local = in{4};
        E = in{5};
      end
    else
      error ('Unknown geometry');
    end
  elseif ismatrix(in{2})
    if size(in{2},1) == 2
      if size(in,2) == 2
        P = in{2};
        N = size(P,2);
      elseif size (in,2) == 3
        P = in{2};
        N = size(P,2);
        C = in{3};
      elseif size (in,2) == 4
        P = in{2};
        N = size(P,2);
        C = in{3};
        local = in {4};
      else
        error ('Illegal number of arguments for matrix');
      end
    else
      error ('Illegal positions');
    end
  else
    error ('Illegal geometry specification');
  end

  % Construct Hamiltonian

  % Compute distances squared
  X = ones(N,1) * P(1,:);
  Y = ones(N,1) * P(2,:);
  % Compute coupling matrix assuming inverse cube law
  H = ((X-X').^2+(Y-Y').^2).^(-3/2);
  H(isinf(H)) = 0;

  % Neighbour couplings only.
  if local ~= 0
    DT = DelaunayTri(P');
    Edges = edges(DT);
    Ne=size(Edges,1);
    DT = DT.Triangulation;
    J = zeros(N);
    for ii=1:Ne
      J(Edges(ii,1),Edges(ii,2)) = H(Edges(ii,1),Edges(ii,2));
    end
    H = J + J';
  end

  % Modify Hamiltonian according to coupling type
  if ~isstr(C)
    error('Quantum spin network coupling type must be a string');
  end
  if strcmpi(C,'XX')
    ;
  elseif strcmpi(C,'H')
    % Heisenberg coupling in first excitation sub-space
    H = H + diag (1/2 * sum(sum(triu(H))) * ones(N,1) - sum(H,2));
  else
    error (['Unknown spin coupling type ' in{2}]);
  end

  % Add engineered network fields
  if ~isnan(E)
    if ~isvector (E)
      error ('Engineered network argument must be a vector');
    end
    if size(E,1) > 1
      E = E';
    end
    if size(E,2) < N
      E = [E ones(1,N-size(E,2))];
    elseif size(E,2) ~= N
      error (sprintf('Engineered ring vector length (%d) must be equal to the number of spins (%d)', size(E,2), N));
    end
    H = H + diag(E);
  end

  % Store data
  obj.N = N;
  obj.H = H;
  obj.pos = P;
  obj.excite = 1;

end