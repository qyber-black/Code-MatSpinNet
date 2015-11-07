function [best,fastest,Results,FailedRuns] = example_bias_hes()
  % Static bias control for information transfer in higher excitation subspaces.
  % in N-ring.

  N   = 7;               % Size
  ring = false;          % Ring?

  %mapping = '1to7_2to6_3to5'; % id string
  %targets = 4;           % Number of targets in M, In, Out
  %M   = { 3, 1, 1, 1 };           % Subspace
  %In  = { '1110000', '1000000', '0100000', '0010000' };   % Input state
  %Out = { '0000111', '0000001', '0000010', '0000100' };   % Output state

  %mapping = 'special'; % id string
  %targets = 3;           % Number of targets in M, In, Out
  %M   = { 2, 1, 1 };           % Subspace
  %In  = { '1100000', '1000000', '0100000' };   % Input state
  %Out = { '0000011', '0000010', '0000100' };   % Output state

  %mapping = 'toffoli'; % id string
  %targets = 7;           % Number of targets in M, In, Out
  %M   = { 1, 1, 2, 1, 2, 3, 3 };           % Subspace
  %In  = { '0010000', '0100000', '0110000', '1000000', '1010000', '1101000', '1110000' };   % Input state
  %Out = { '0000001', '0000010', '0000011', '0000100', '0000101', '0000111', '0001110' };   % Output state

  N   = 11;               % Size
  mapping = 'toffoli_rev'; % id string
  targets = 7;           % Number of targets in M, In, Out
  M   = { 1, 1, 2, 1, 2, 3, 3 };           % Subspace
  In  = { '00100000000', '01000000000', '01100000000', '10000000000', '10100000000', '11010000000', '11100000000' };   % Input state
  Out = { '00000000100', '00000000010', '00000000110', '00000000001', '00000000101', '00000000111', '00000001011' };   % Output state

  C_symmetry = false;  % Symmetric controls around central spin?

  maxT = 199;           % Max init time - 1
  maxB = 1000;         % Max init bias

  diagonals = false;      % simple controls on the diagonal (only with single target!)
  trust_region = false; % trust-region or quasi newton?

  min_err = 0.01;      % Largest error acceptable for shortest solution.
  repeats = 100;       % How many restarts?

  if C_symmetry
    if ring
      fname = sprintf('ring_%d_sym_%s',N,mapping);
    else
      fname = sprintf('chain_%d_sym_%s',N,mapping);
    end
  else
    if ring
      fname = sprintf('ring_%d_%s',N,mapping);
    else
      fname = sprintf('chain_%d_%s',N,mapping);
    end
  end

  % States
  Bin = dec2base([0:2^N-1],2);
  for t =1:targets
    In{t} = fliplr(In{t});
    Out{t} = fliplr(Out{t});
    Ind = find(sum(Bin=='1',2)==M{t});
    for k = 1:length(Ind)
      if strcmp(Bin(Ind(k),:),In{t}) ~= 0
      IN{t} = k;
      end
      if strcmp(Bin(Ind(k),:),Out{t}) ~= 0
        OUT{t} = k;
      end
    end
    disp(sprintf('%d: From %s (%d) to %s (%d)', M{t}, fliplr(Bin(Ind(IN{t}),:)), IN{t}, fliplr(Bin(Ind(OUT{t}),:)), OUT{t}));
  end

  % XX-ring in higher excitation subspaces and Z controls
  JX = diag(ones(1,N-1),1) + diag(ones(1,N-1),-1);
  if ring
    JX(1,N) = 1;
    JX(N,1) = 1;
  end
  JZ = zeros(N,N);
  [H0,C] = HES(JX,JX,JZ,M);
  N_C = size(C,2);

  % Optimisation loop
  figure(1);
  clf;
  drawnow ();
  disp('Solving...');
  best = 0;
  fastest = 0;
  Results = cell(1,repeats);
  Err = zeros(1,repeats);
  FailedRuns = {};

  for rep = 1:repeats
    failed = true;

    % Repeats over failed optimisations
    while failed
      % Initial value
      maxBB = maxB * rand(1);
      x0 = rand (1,N_C);
      x0(1,1:N_C) = (x0(1,1:N_C) / max(x0(1,1:N_C))) * maxBB;
      x0(1,N_C+1) = 1 + maxT * rand(1);
      % Optimisation
      failed = false;
      try
        options = optimoptions('fminunc');
        if trust_region
          options.Algorithm = 'trust-region';
          options.CheckGradients = false;
          options.MaxFunctionEvaluations = 1000 * N_C;
          options.MaxIterations = 1000 * N_C;
          options.SpecifyObjectiveGradient = true;
          options.SubproblemAlgorithm = 'factorization';
        else
          options.Algorithm = 'quasi-newton';
          options.MaxFunctionEvaluations = 1000 * N_C;
          options.MaxIterations = 1000 * N_C;
          options.SubproblemAlgorithm = 'factorization';
        end
        [Results{rep}.x,Results{rep}.err,Results{rep}.exit_flag,Results{rep}.output] = fminunc (@eval_static, x0, options);
        Results{rep}.err = eval_static(Results{rep}.x); % Ensure fidelity error is reported
        Results{rep}.x0 = x0;
      catch ME
        display('Failed optimisation run with x0 = ');
        disp(x0);
        getReport(ME)
        failed = true;
        fidx = size(FailedRuns,2) + 1;
        FailedRuns{fidx}.x0 = x0;
        FailedRuns{fidx}.reason = getReport(ME);
      end
      if ~failed
        if Results{rep}.exit_flag < 0
          display(sprintf('Failed exit_flag = %d run with x0 = ', Results{rep}.exit_flag));
          disp(x0);
          fidx = size(FailedRuns,2) + 1;
          FailedRuns{fidx}.x0 = x0;
          FailedRuns{fidx}.reason = sprintf('Exist flag %d',Results{rep}.exit_flag);
          failed = true;
        elseif Results{rep}.err < 0
          display('Negative error run with x0 = ');
          disp(x0);
          fidx = size(FailedRuns,2) + 1;
          FailedRuns{fidx}.x0 = x0;
          FailedRuns{fidx}.reason = sprintf('Error %g negative',Results{rep}.err);
          failed = true;
        elseif Results{rep}.err > 0.5
          display('Too large error run with x0 = ');
          disp(x0);
          fidx = size(FailedRuns,2) + 1;
          FailedRuns{fidx}.x0 = x0;
          FailedRuns{fidx}.reason = sprintf('Error %g too large',Results{rep}.err);
          failed = true;
        end
      end
    end

    % Plot
    plot_result(rep, 2,4,1,5,'Current solution');

    % Best solution
    if best == 0
      best = rep;
      plot_result(best, 2,4,2,6,'Best solution');
    else
      if Results{rep}.err < Results{best}.err
        best = rep;
        plot_result(best, 2,4,2,6,'Best solution');
      end
    end

    % Fastest solution
    if Results{rep}.err < min_err
      if fastest == 0
        fastest = rep;
        plot_result(fastest,2,4,3,7,'Fastest solution');
      elseif Results{rep}.x(end) < Results{fastest}.x(end)
        fastest = rep;
        plot_result(fastest,2,4,3,7,'Fastest solution');
      end
    end

    % Error vs time
    subplot(2,4,4);
    hold on;
    plot(log10(Results{rep}.err),Results{rep}.x(end),'*b');
    hold off;
    xlabel('log10(Error)');
    ylabel('Time');
    title('log10(Error) vs. Time');
    axis tight;

    % Error Histogram
    subplot(2,4,8);
    Err(rep) = log10(Results{rep}.err);
    histogram(Err(1:rep),'FaceColor',[0 0 1],'FaceAlpha',1);
    title(sprintf('log10(Error) Histogram, %d runs', rep));
    axis tight;

    drawnow();
    refresh();

  end

  % Final plot
  clf
  % Best solution
  if best > 0
    plot_result(best, 2,3,1,4,'Best solution');
  end
  % Fastest solution
  if fastest > 0
    % Setup trace plots
    plot_result(fastest, 2,3,2,5,'Fastest solution');
  end

  % Error vs time
  subplot(2,3,3);
  hold on;
  Times = arrayfun(@(x) (Results{x}.x(end)),[1:size(Results,2)]);
  Err = log10(arrayfun(@(x) (Results{x}.err),[1:size(Results,2)]));
  plot(Err, Times, '*b');
  xlabel('log(Error)');
  ylabel('Time');
  title('log(Error) vs. Time');
  axis tight;

  % Error Histogram
  subplot(2,3,6);
  histogram(Err,'FaceColor',[0 0 1],'FaceAlpha',1);
  title(sprintf('log(Error) Histogram over %d runs', size(Results,2)));
  axis tight;

  display(sprintf('Failed runs: %d\n', size(FailedRuns,2)));

  drawnow();
  refresh();

  savefig(fname);
  save(fname,'best','fastest','Results','FailedRuns');

  % Plot results
  function plot_result (run, X, Y, trace_fig, bias_fig, str)
    % Plot traces
    subplot(X,Y,trace_fig);
    plot_time = [0:Results{run}.x(end)/1000:Results{run}.x(end)];
    plot_ctrl = [];
    plot_nat = [];
    for ta = 1:targets
      H = H0{ta};
      for l = 1:N_C
        H = H + Results{run}.x(l) * C{ta,l};
      end
      for t = 1:size(plot_time,2)
        U = expm(-1i*H*plot_time(t));
        plot_ctrl(t) = abs(U(OUT{ta},IN{ta}))^2;
        U = expm(-1i*H0{ta}*plot_time(t));
        plot_nat(t) = abs(U(OUT{ta},IN{ta}))^2;
      end
      plot(plot_time, plot_ctrl); %, plot_time, plot_nat);
      axis([0 abs(plot_time(end)) 0 1]);
      hold on;
      %plot(plot_time(end),plot_ctrl(end),'*r');
      %plot(plot_time(end),plot_nat(end),'*g');
    end
    hold off;
    xlabel('time in 1/J')
    ylabel(sprintf('p_{%i%i}',IN{1},OUT{1}));
    title(sprintf('%s, N=%d, e=%.6g', str, N, Results{run}.err));
    % Plot bias
    subplot(X,Y,bias_fig);
    if C_symmetry
      if N/2 == floor(N/2)
        bar(1:N,[Results{run}.x(1:N_C) Results{run}.x(N_C:-1:1)],'b');
      else
        bar(1:N,[Results{run}.x(1:N_C) Results{run}.x(N_C-1:-1:1)],'b');
      end
    else
      bar(1:N_C,Results{run}.x(1:N_C),'b');
    end
    hold on;
    if C_symmetry
      if N/2 == floor(N/2)
        bar(1:N,[Results{run}.x0(1:N_C) Results{run}.x0(N_C:-1:1)],0.2,'g');
      else
        bar(1:N,[Results{run}.x0(1:N_C) Results{run}.x0(N_C-1:-1:1)],0.2,'g');
      end
    else
      bar(1:N_C,Results{run}.x0(1:N_C),0.2,'g');
    end
    hold off;
    xlabel('Spin #');
    ylabel('Bias');
    title(sprintf('Bias |%d>-|%d>, T=%.6g, T0=%.6g',IN{1},OUT{1},Results{run}.x(end),Results{run}.x0(end)));
    axis tight;
  end

  function [err,grad] = eval_static (x)
    if nargout > 1
      [err,grad] = eval_static_single(x,1);
      for t = 2:targets
        [e,g] = eval_static_single(x,t);
        err = err + e;
        grad = grad + g;
      end
    else
      err = eval_static_single(x,1);
      for t = 2:targets
        err = err + eval_static_single(x,t);
      end
    end
  end

  function [err,grad] = eval_static_single (x,t)
    T = abs(x(end)); % ensure T >=0
    H = H0{t};
    for l = 1:N_C
      H = H + x(1,l) * C{t,l};
    end
    %U = expm(-1i*T*H);
    [V,D] = eig(-1i*T*H);
    U = V * diag(exp(diag(D))) / V;
    phi = U(OUT{t},IN{t});
    err = 1-abs(phi)^2; % FIXME: minimise time?
    if nargout > 1,
      grad = zeros (N_C,1); % derivative of the error wrt f(mk)
      for m=1:N_C
        [~,dU]  = dexpma(-1i*T*H,-1i*T*C{t,m});
        grad(m) = -2*real(dU(OUT{t},IN{t})*conj(phi));
      end
      % dF/dT = -2 Re <p2|-iH*U|p1> conj(<p2|U|p1>)
      grad(N_C+1) = -2*imag(H(OUT{t},:)*U(:,IN{t})*conj(phi));
    end
  end

 function [err,grad] = eval_static_dephasing(x,t,c)
    % c are the eigenvalues of the dephasing operator
    T = abs(x(end)); % ensure T >=0
    H = H0{t};
    for l = 1:N_C
      H = H + x(1,l) * C{t,l};
    end
    %U = expm(-1i*T*H);
    [V,D] = eig(H);  d = diag(D); 
    iomega = 1i*(d*ones(1,length(d))-ones(length(d),1)*d');
    gamma  = 1/2*(c*ones(1,length(c))-ones(length(c),1)*d').^2;
    
    rho0 = zeros(size(H)); rho0(IN{t},IN{t})   = 1;
    rho1 = zeros(size(H)); rho1(OUT{t},OUT{t}) = 1;
    err  = 1-trace(rho1*(-T*(iomega+gamma).*(V'*rho0*V)));
    
  end

  function [F,dF] = dexpma(M,dM,v)
    % [F,dF] = dexpma(M,dM,v)
    %
    % Gradient helper
    N   = length(M);
    AM  = [M zeros(N);dM M];
    if exist('v','var')
      PSI = expv(AM,[v;v]);
      F   = PSI(1:N);
      dF  = PSI(N+1:end);
    else
      PSI = expm(AM);
      F   = PSI(1:N,1:N);
      dF  = PSI(N+1:2*N,1:N);
    end
  end

  function [Hs,Cs] = HES(Jx,Jy,Jz,SubSpaces)
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

    % If H is the full 2^N dimensional Hamiltonian for the N spin system
    % I is the index vector so that H(I,I) is block-diagonal with the kth
    % block corresponding to the (k-1)th excitation subspace
    ind = [];
    if diagonals
      for s = 1:length(SubSpaces)
        ind = ESubspace(N,SubSpaces{s});
        Hs{s} = H(ind,ind);
        % Controls
        K = size(ind,2)
        for m = 1:K
          Cs{s,m} = zeros(K,K);
          Cs{s,m}(m,m) = 1;
          Cs{s,m}
        end
      end
    else
      for s = 1:length(SubSpaces)
        ind = ESubspace(N,SubSpaces{s});
        Hs{s} = H(ind,ind);
        % Controls
        for m = 1:N
          K = kron(kron(eye(2^(m-1)),Z),eye(2^(N-m)));
          Cs{s,m} = K(ind,ind);
          Cs{s,m} = (diag(ones(size(Cs{s,m},1),1))-Cs{s,m})/2;
          % diag(Cs{s,m})'
        end
        if C_symmetry
          for m = 1:floor(N/2)
            Cs{s,m} = Cs{s,m} + Cs{s,N+1-m};
          end
        end
      end
      if C_symmetry
        for s = 1:length(SubSpaces)
          for m=1:ceil(N/2)
            CC{s,m} = Cs{s,m};
          end
        end
        clear Cs;
        Cs = CC;
      end
    end
  end

  function ind = ESubspace(N,k)
    % computes indices I of k excitation subspace for spin chain of length N
    % If H = GenerateNspinH(...), H(I,I) = k-excitation subspace Hamiltonian
    ind = find(sum(dec2base([0:2^N-1],2)=='1',2)==k);
    ind = fliplr(ind');
  end

end
