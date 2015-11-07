function create_hes (obj, in)
  % create_hes (obj, in) - Create a higher excitation subspace spin network
  %
  % Construct a quantum spin network considering the selected excitation
  % sub-spaces.
  %
  %   in     - Cell array describing 2D network as follows:
  %     in{1} - if 1 matrix assume isotropic Heisenberg
  %             if 3 matrices in cell then anisotripic Heisenberg
  %             if 2 matrices in cell assume XY interaction
  %             if 1 matrices in cell isotropic Heisenberg
  %     in{2} - vector of excitation subspaces to represent,
  %             in which order (0:k - from 1 to k <=N)
  %     in{3} - row vector of biases on the spins (optional)

  % Check arguments
  if size(in,2) < 2 || size(in,2) > 3
    error ('Wrong number of parameters for higher-excitation spin network');
  end

  % Coupling matrices
  J = in{1};
  if ~isa(J,'cell'),      %%% if single matrix assume isotropic Heisenberg
    Jx = J;
    Jy = J;
    Jz = J;
  elseif (length(J)==3),  %%% if 3 matrices then anisotripic Heisenberg
    Jx = J{1};
    Jy = J{2};
    Jz = J{3};
  elseif (length(J)==2),  %%% if 3 matrices assume XY interaction
    Jx = J{1};
    Jy = J{2};
    Jz = zeros(size(J{1}));
  elseif (length(J)==1),  %%% assume isotropic Heisenberg
    Jx = J{1};
    Jy = Jx;
    Jz = Jx;
  else
    error('cannot interpret coupling matrices')
  end

  % Size of network
  N = size(Jx,1);

  % Excitation subspaces
  es = in{2};

  % Construct positions: just on a chain
  P = [ ([0:N-1] / (N-1) * 2) - 1; zeros(1,N) ];

  % Construct Hamiltonian:
  % constructs the full 2^N x 2^N Hamiltonian for an array of N spins with
  % anisotropic Heisenberg coupling between arbitary spins given by the
  % coupling J={Jx,Jy,Jz} where Jx, Jy, Jz are symmetric coupling matrices

  % FIXME: alternative calculation of H for sub-spaces only from graph paper

  % define pauli matrices
  X = [0 1;1 0];
  Y = [0 -i;i 0];
  Z = [1 0;0 -1]; %%Z=[-1 0;0 1];
  % generate N-qubit Hamiltonian
  H=0;
  [M,N]=size(Jx);
  for m=1:M
    for n=m+1:N
      if Jx(m,n) ~=0
        H=H+Jx(m,n)*kron(kron(kron(eye(2^(m-1)),X),eye(2^(n-m-1))),kron(X,eye(2^(N-n))));
      end
      if Jy(m,n) ~=0
        H=H+Jy(m,n)*kron(kron(kron(eye(2^(m-1)),Y),eye(2^(n-m-1))),kron(Y,eye(2^(N-n))));
      end
      if Jz(m,n) ~=0
        H=H+Jz(m,n)*kron(kron(kron(eye(2^(m-1)),Z),eye(2^(n-m-1))),kron(Z,eye(2^(N-n))));
      end
    end
  end
  H = H/2;

  if size(in,2) == 3
    b = in{3};
    if size(b,2) ~= N
      error('Bias vector does not match number of spins')
    end
    for m = 1:N
      H = H + b(1,m) * kron(kron(eye(2^(m-1)),Z),eye(2^(N-m)));
    end
  end

  % If H is the full 2^N dimensional Hamiltonian for the N spin system
  % I is the index vector so that H(I,I) is block-diagonal with the kth
  % block corresponding to the (k-1)th excitation subspace
  I = [];
  for s = es
    I = [I,ESubspace(N,s)];
  end

  % Store data
  obj.N = N;
  obj.H = H(I,I);
  obj.pos = P;
  obj.excite = max(es);

  function I = ESubspace(N,k)
    % computes indices I of k excitation subspace for spin chain of length N
    % If H = GenerateNspinH(...), H(I,I) = k-excitation subspace Hamiltonian
    I = find(sum(dec2base([0:2^N-1],2)=='1',2)==k);
    I = fliplr(I');
  end
  % matlab-only version
  % I=find(sum(dec2base([0:2^N-1],2),2)-48*N==k);
  %%% -48*N is necessary to account for the fact that dec2base(X,2) produces
  %%% a char array of length N with entries 0 and 1.  Since the ascii values
  %%% of 0 and 1 are 48 and 49, respectively, sum(dec2base(X,2),2)-48*N counts
  %%% the number of 1's in the binary representation of X, assuming it is a
  %%% string of length N.

end
