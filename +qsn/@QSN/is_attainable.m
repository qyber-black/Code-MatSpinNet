function [p,q,theta,OK,err,info] = is_attainable (obj,ispin,jspin,x,varargin)
  % obj . is_attainable (ispin, jspin, x, shortest) - Check if probability bound is attainable (rings)
  %
  % This function checks whether, in a spin network ring with XX or Heisenberg
  % couplings, p_{max}(ispin,jspin) can be achieved, i.e. whether the probability
  % of transfer of the excitation from ispin to jspin can reach its upper bound.
  %
  % For this we use the simultaneous diophantine approximation.
  %
  % If shortest is given and not 0, then find shortest vector for simultaneous
  % Diophantine approximation.
  %
  % This latter diophantine approximation utilizes x as scaling parameters.
  % If x is a single number, then uniform scaling will be used. If x is a
  % (d+1)-dimensional vector [y s] then y is the non-uniform scaling vector and
  % s is a scalar scaling parameter. The dimension d=ceil((N+1)/2)-2 must be
  % entered manually as the number of components of x in this case.
  %
  % On return, we have the near optimal simultaneous diophantine
  % approximation \theta q ~ p, where \theta is constructed from the
  % translation on the torus.  The accuracy of the approximation as x goes to
  % infinity and s goes to 0 can be checked by computing the error as
  %
  % error = theta*q - p
  %
  % The above error is indeed the error relevant to accuracy of the attainablity.
  % The above is denoted as \epsilon_{\mathrm{Da}} in the paper.
  %
  % The even/odd specifications on the vector p of numerators are in even_odd,
  % a vector of the same length as p made up of 0's or \pm 1's depending whether
  % the corresponding component of p should be even/odd, resp.
  %
  % All relevant results to see "where we stand" are consolidated in the
  % info, also on output, with the p in the first column, whether the p should
  % be even (==0) or odd (==1) in the second column, whether p meets the
  % even/odd condition in the third column (1=OK, 0=not), and the corresponding
  % errors theta * q - p in the fourht column as the rounded value of the log
  % in base 10.
  %
  % Output err is the error theta*q-p and OK indicates if the even/odd conditions
  % are met by p.
  %

  if size(varargin,2) < 1
    shortest = 0;
  else
    shortest = varargin{1};
  end

  if strcmpi (obj.type, 'ring')
    test_ring = qsn.QSN('ring',obj.N,'XX');
    d = norm(obj.H - test_ring.H,1);
  else
    d = 0;
  end
  if d < obj.epsilon()
    % For unbiased ring
    %
    Ntilde = ceil ((obj.N + 1) / 2);
    Nbar = Ntilde - 2;
    
    % We enter \theta_{k,\ell} from the translation on the torus,
    % assuming \ell=k+1:
    %
    theta = transpose (sin((pi*(2*[1:Nbar]+1))/obj.N) * 2 / sin(pi/obj.N));
    %
    % We derive the eigenstructure of the Hamiltonian matrix;
    % 
    % eigenvalues
    %
    e = 2*cos(2*pi*[0:obj.N-1]/obj.N);
    %
    % eigenvectors
    %
    w = cellfun (@(k) transpose (exp(2*pi*i*k*[0:obj.N-1]/obj.N) / sqrt(obj.N)), num2cell(1:obj.N), 'UniformOutput', false);
    %
    % projection operators  
    %
    W = cellfun (@(w) w * w', w, 'UniformOutput', false);
    if mod(obj.N,2)
        Pi = arrayfun (@(k) W{k} + W{obj.N-k}, (1:(obj.N-1)/2), 'UniformOutput', false);
        Pi = {W{end}, Pi{1:end}};
        e = e(1:(obj.N+1)/2);
    else
        Pi = arrayfun (@(k) W{k} + W{obj.N+1-k}, (1:obj.N/2), 'UniformOutput', false);
        e = e(1:obj.N/2);
    end
    % check decomposition is correct
    H = zeros(obj.N); 
    for k=1:length(Pi)
        H = H + e(k)*Pi{k}; 
    end
    if norm(H-obj.H)>1e-12,
        error('we have a problem')
    end
    % Signs
    sg = real(cell2mat (arrayfun (@(k) sign(Pi{k}(ispin,jspin)), (1:ceil(obj.N/2)), 'UniformOutput', false)));
    %
    % Finally the simultaneous diophantine aproximation of \theta:
    %
    [p,q,short] = sim_dio_approx (theta, x, shortest);
    keyboard
  else
    % General

    % Eigenvectors and unique eigenvalues
    [V,E] = eig(obj.H);
    E = diag(E);
    % For consistency with ring reverse order
    %E = E(size(E,1):-1:1);
    %V = V(:,size(V,2):-1:1);
    [~,indE,induE] = unique(floor(E/obj.epsilon()));
    uE = E(indE);
    % -> uE = E(indE); E = uE(induE);
    Ntilde = size(uE,1);

    % Projectors onto eigenspaces
    VV = cellfun(@(k) V(:,k) * V(:,k)', num2cell(1:size(V,2)), 'UniformOutput', false);
    Pi = cell(1,Ntilde);
    for k = 1:size(VV,2)
      ind = induE(k);
      if isempty (Pi{ind})
        Pi{ind} = VV{k};
      else
        Pi{ind} = Pi{ind} + VV{k};
      end
    end

    % Signs
    sg = cell2mat (cellfun (@(k) sign(Pi{k}(ispin,jspin)), num2cell(1:Ntilde), 'UniformOutput', false));

    % (Non-)Dark states
    indNDS = find(sg ~= 0);
    sg = sg(indNDS);
    uE = uE(indNDS);

    % omega_mn - two equal signs
    M = NaN;
    N = NaN;
    for m = size(uE,1):-1:2
      for n = m-1:-1:1
        if sg(m) == sg(n)
          M = m;
          N = n;
          break
        end
      end
      if ~isnan(M)
        break
      end
    end
    if isnan(M)
      error('Unecpected result: omega_mn not found');
    end

    % Get omega_mn and resort rest of eigenvalues to calcuate remainign omega_lk
    omega_mn = (uE(M) - uE(N)) / pi;
    uE = uE([1:m-1 m+1:size(uE,1)]);
    sg = sg([1:m-1 m+1:size(uE,1)]);

    % We enter \theta_{k,\ell} from the translation on the torus,
    % assuming \ell=k+1:
    %
    theta = 2*(diff(uE) / pi) / omega_mn;

    % Next, the simultaneous diophantine aproximation of \theta:
    %
    [p,q,short] = sim_dio_approx (theta, x, shortest);
  end

  % Now we define the vector with 0 or \pm 1 components,
  % indicating that the numerators must even or odd, resp., for attainability.
  % We call the vector
  %
  % (1/2)(s_k-s_\ell)
  %
  % the vector of specifications on p, specp.
  %
  even_odd = abs (transpose (diff(sg) / 2)); % 0 = even; 1 = odd;
  %
  % Check results
  %
  match = (even_odd == abs(mod(p,2)));
  OK = (sum(match) == size(theta,1));
  err = theta*q - p;
  info = [p even_odd match round(log(abs(err)))];
end