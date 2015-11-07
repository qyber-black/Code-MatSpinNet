function [Results,best,fastest,FailedRuns,info] = bias_control(obj,in,out,T,B,readout,min_err,repeats,symm,initT,bias_init,noise,v)
  % [Results,best,fastest] = bias_control(obj,in,out,T,B,readout,min_err,repeats,symm,initT,bias_init,noise)
  % Find static bias control for QSN
  %
  % in         input spin (row vector for more than one)
  % out        readout spin (row vector for more than one)
  % T          Fixed time or NaN to find time
  % B          Maximum bias or NaN for no constraint
  % readout    Readout window (or empty to ignore)
  %                     readout(1) Time window for readout (0 - do not consider window)
  %                     readout(2) Perturbation (or 0 for none, time window required)
  %                     readout(3) Samples for expectation over perturbation (if perturbation not 0)
  % min_err    Largest acceptable error for fastest solution
  % repeats    How many repeats?
  % symm       enforce symmetric controls for chain, ring, xring only
  % initT      empty or initT(1) - maximum time to search for peaks
  %                     initT(2) - number of samples between 0 and initT(1) to search for peaks
  %                              - any further entries are added to the initial times selected
  % bias_init  initialise biases with peak/trough/constant instead of random (0)?
  % noise      Add level of noise to initial values
  % v          Verbose?

  % Verbose by default
  if ~exist('v','var')
    v = 1;
  end

  % Result figure
  if v ~= 0
    figure(1);
    colormap('default');
  end

  rng('shuffle');

  info.args.obj = obj;
  info.args.in = in;
  info.args.out = out;
  info.args.T = T;
  info.args.B = B;
  info.args.readout = readout;
  info.args.min_err = min_err;
  info.args.repeats = repeats;
  info.args.symm = symm;
  info.args.initT = initT;
  info.args.bias_init = bias_init;
  info.args.noise = noise;

  if size(in,1) ~= size(out,1)
    error('input spin vector must match output spin vector')
  end

  % Setup symmetric controls
  if symm == 0
    Cspin = 1:obj.N;
  else
    if size(in,1) > 1
      error('symmetry currently only implemented for single excitation sub-space')
    end
    if strcmpi(obj.type,'chain') || strcmpi(obj.type,'ring')
      if in < out
        f = in;
        l = out;
      else
        f = out;
        l = in;
      end
      gap = floor((l-f)/2);
      for k = 0:gap
        Cspin(f+k) = k+1;
        Cspin(l-k) = k+1;
      end
      for k = 1:ceil((obj.N-(l-f)-1)/2)
        if f-k < 1
          Cspin(obj.N+(f-k)) = gap+k+1;
        else
          Cspin(f-k) = gap+k+1;
        end
        if l+k > obj.N
          Cspin(l+k-obj.N) = gap+k+1;
        else
          Cspin(l+k) = gap+k+1;
        end
      end
    elseif strcmpi(obj.type,'xring')
      Nring = max(find(sum(obj.H) == 3));
      inring = in;
      while inring > Nring
        inring = min(find(obj.H(inring,:) == 1));
      end
      outring = out;
      while outring > Nring
        outring = min(find(obj.H(outring,:) == 1));
      end
      if inring < outring
        f = inring;
        l = outring;
      else
        f = outring;
        l = inring;
      end
      gap = floor((l-f)/2);
      for k = 0:gap
        Cspin(f+k) = k+1;
        Cspin(l-k) = k+1;
      end
      for k = 1:ceil((Nring-(l-f)-1)/2)
        if f-k < 1
          Cspin(Nring+(f-k)) = gap+k+1;
        else
          Cspin(f-k) = gap+k+1;
        end
        if l+k > Nring
          Cspin(l+k-Nring) = gap+k+1;
        else
          Cspin(l+k) = gap+k+1;
        end
      end
      Cpos = max(Cspin)+1;
      for k = 1:Nring
        fs = max(find(obj.H(k,:) == 1));
        if fs > Nring
          ls = k;
          ns = fs;
          while ns > ls
            ls = ns;
            ns = max(find(obj.H(ns,:) == 1));
          end
          gap = floor((ls-fs)/2);
          for k = 0:gap
            Cspin(fs+k) = Cpos;
            Cspin(ls-k) = Cpos;
            Cpos = Cpos + 1;
          end
        end
      end
    else
      error('Symmetry not implemented');
    end
  end

  % Initialise symmetric ranges in network with peaks/troughs/constant?
  Ranges = [];
  if bias_init ~= 0
    if size(in,1) > 1
      error('symmetric initialisation currently only implemented for single excitation sub-space')
    end
    if strcmpi(obj.type,'chain') || strcmpi(obj.type,'ring')
      if in < out
        f = in;
        l = out;
      else
        f = out;
        l = in;
      end
      Ranges = [Ranges; f l];
      l = l + 1;
      if l > obj.N
        l = 1;
      end
      f = f - 1;
      if f < 1
        f = obj.N;
      end
      Ranges = [Ranges; l f];
    elseif strcmpi(obj.type,'xring')
      Nring = max(find(sum(obj.H) == 3));
      inring = in;
      while inring > Nring
        inring = min(find(obj.H(inring,:) == 1));
      end
      outring = out;
      while outring > Nring
        outring = min(find(obj.H(outring,:) == 1));
      end
      if inring < outring
        f = inring;
        l = outring;
      else
        f = outring;
        l = inring;
      end
      Ranges = [Ranges; f l];
      l = l + 1;
      if l > obj.N
        l = 1;
      end
      f = f - 1;
      if f < 1
        f = obj.N;
      end
      Ranges = [Ranges; l f];
      for k = 1:Nring
        fs = max(find(obj.H(k,:) == 1));
        if fs > Nring
          ls = k;
          ns = fs;
          while ns > ls
            ls = ns;
            ns = max(find(obj.H(ns,:) == 1));
          end
          Ranges = [Ranges; fs ls];
        end
      end
    else
      error('bias_init not implemented');
    end
  end

  % Spin l is controlled by bias in C{CSpin(l)}
  C = cell(1,max(Cspin));
  for l = 1:size(C,2)
    C{l} = zeros(obj.N,obj.N);
  end
  for l = 1:size(Cspin,2)
    C{Cspin(l)}(l,l) = 1;
  end

  info.controls = C;

  % For testing
  %Cspin
  %size(Cspin)
  %CC = C{1};
  %size(CC)
  %for l = 2:size(C,2)
  %  size(C{l})
  %  CC = CC + l * C{l};
  %end
  %heatmaptext(CC);
  %pause

  init_times = [];
  if isnan(T)
    if size(initT,2) > 1
      if size(in,1) > 1
        error('Time from chian currenlty only implemented for single-excitation sub-spaces');
      end
      if strcmpi(obj.type,'chain')
        dist = abs(in-out)+1;
      elseif strcmpi(obj.type,'ring')
        dist = [abs(in-out)+1, obj.N-abs(in-out)+2];
      elseif strcmpi(obj.type,'xring')
        dist = 0;
        Nring = max(find(sum(obj.H) == 3));
        inring = in;
        while inring > Nring
          dist = dist + 1;
          inring = min(find(obj.H(inring,:) == 1));
        end
        outring = out;
        while outring > Nring
          dist = dist + 1;
          outring = min(find(obj.H(outring,:) == 1));
        end
        if inring == outring
          dist = abs(in - out) + 1;
        else
          dist = dist + [abs(inring - outring)+1, obj.N-abs(inring-outring)+2];
        end
      else
        error('Distance for chain time not implemented');
      end
      % Find times with highest probabilities for related chains
      % FIXME: use time estimates here!
      cT = 0:initT(1)/(initT(2)-1):initT(1);
      for d = dist
        chain = qsn.QSN('chain',  d);
        P = chain.prob(cT);
        pit = arrayfun(@(n) P{n}(1,chain.N), 1:length(P));
        [~,ti] = findpeaks(pit, 'MinPeakHeight', .5, 'SortStr', 'descend');
        init_times = [init_times, cT(ti)];
      end
      if size(initT,2) > 2
        init_times = [init_times, initT(3:end) ];
      end
      init_times = sort(init_times);
      time_index = 1;
    end
  end

  % Setup control problem
  if size(readout,2) ~= 0
    if readout(1) == 0
      readout = [];
    end
    if size(readout,2) > 1
      if readout(2) == 0
        readout = readout(1);
      elseif size(readout,2) ~= 3
        error('illegal readout parameter');
      end
    end
  end
  if size(in,1) == 1
    IN = in;
    OUT = out;
  else
    % Find states corresponding to HES transfer
    error('not implemented');
  end
  if size(readout,2) == 0
    ctrl = qsn.Control (obj, 'state_transfer', C,  IN, OUT, T, B);
  else
    if size(readout,2) == 1
      ctrl = qsn.Control (obj, 'state_transfer', C, IN, OUT, T, B, readout(1));
    else
      ctrl = qsn.Control (obj, 'state_transfer', C, IN, OUT, T, B, readout(1), readout(2), readout(3));
    end
  end

  % Optimisation loop
  if v ~= 0
    ctrl
    disp('Solving...');
    clf;
  end
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
      if isnan(B)  % Strength
        maxBB = 1000 * rand(1);
      else
        maxBB = B;
      end
      lC = length(C);
      if bias_init == 0
        x0 = rand (1,lC);
      else
        % Initialise Ranges with peaks/troughs/constants
        sw = rand(1,size(Ranges,1));
        x0 = nan(1,lC);
        for r = 1:size(Ranges,1)
          f = Ranges(r,1);
          l = Ranges(r,2);
          if f > l
            gap = obj.N-l + f + 1;
          else
            gap = l-f + 1;
          end
          if sw(1,r) < 1/3
            % Peak between f and l
            val = 1 - abs([1:gap] - (gap+1)/2) - (1-(gap+1)/2);
          elseif sw(1,r) < 2/3
            % Constants between f and l
            val = ones(1,gap);
          else
            % Trough between f and l
            val = floor(1+abs([1:gap] - (gap+1)/2));
          end
          if f > l
            x0(1,Cspin([1:f l:obj.N])) = val;
          else
            x0(1,Cspin(f:l)) = val;
          end
        end
        ind = find(isnan(x0) == 1);
        if ~isempty(ind)
          x0(ind) = rand(1,size(ind,2));
        end
      end
      x0(1,1:lC) = (x0(1,1:lC) / max(x0(1,1:lC)) + rand(1,lC) * noise) * maxBB;
      if isnan(T) % Time
        if size(init_times,2) == 0
          if size(readout,2) == 0
            x0(1,lC+1) = 1 + 99 * rand(1);
          else
            x0(1,lC+1) = 1 + (readout(1)+1) * 99 * rand(1);
          end
        else
          x0(1,lC+1) = init_times(time_index);
          time_index = time_index + 1;
          if time_index > size(init_times,2)
            time_index = 1;
          end
        end
      end

      % Optimisation
      failed = false;
      try
        if isnan(T) && size(readout,2) == 0
          [x,Results{rep}.err,Results{rep}.exec_time,Results{rep}.exit_flag,Results{rep}.output,x0] = ctrl.solve_static (qsn.Opt ('fmingrad'), x0); % FIXME: make Opt selectable
        else
          [x,Results{rep}.err,Results{rep}.exec_time,Results{rep}.exit_flag,Results{rep}.output,x0] = ctrl.solve_static (qsn.Opt ('fmin'), x0); % FIXME: make Opt selectable
        end
      catch ME
        if v ~= 0
          display('Failed optimisation run with x0 = ');
          disp(x0);
          getReport(ME)
        end
        failed = true;
        fidx = size(FailedRuns,2) + 1;
        init_mat = zeros(obj.N,obj.N);
        for num = 1:size(C,2)
          init_mat = init_mat + x0(1,num) * C{num};
        end
        FailedRuns{fidx}.init_bias = diag(init_mat);
        if isnan(T)
          FailedRuns{fidx}.init_time = x0(1,lC+1);
        else
          FailedRuns{fidx}.init_time = T;
        end
        FailedRuns{fidx}.reason = getReport(ME);
      end
      if ~failed
        if Results{rep}.exit_flag < 0
          if v ~= 0
            display(sprintf('Failed exit_flag = %d run with x0 = ', Results{rep}.exit_flag));
            disp(x0);
          end
          fidx = size(FailedRuns,2) + 1;
          init_mat = zeros(obj.N,obj.N);
          for num = 1:size(C,2)
            init_mat = init_mat + x0(1,num) * C{num};
          end
          FailedRuns{fidx}.init_bias = diag(init_mat);
          if isnan(T)
            FailedRuns{fidx}.init_time = x0(1,lC+1);
          else
            FailedRuns{fidx}.init_time = T;
          end
          FailedRuns{fidx}.reason = sprintf('Exist flag %d',Results{rep}.exit_flag);
          failed = true;
        elseif Results{rep}.err < 0
          if v ~= 0
            display('Failed error run with x0 = ');
            disp(x0);
          end
          fidx = size(FailedRuns,2) + 1;
          init_mat = zeros(obj.N,obj.N);
          for num = 1:size(C,2)
            init_mat = init_mat + x0(1,num) * C{num};
          end
          FailedRuns{fidx}.init_bias = diag(init_mat);
          if isnan(T)
            FailedRuns{fidx}.init_time = x0(1,lC+1);
          else
            FailedRuns{fidx}.init_time = T;
          end
          FailedRuns{fidx}.reason = sprintf('Error %g negative',Results{rep}.err);
          failed = true;
        end
      end
    end

    bias_mat = zeros(obj.N,obj.N);
    init_mat = zeros(obj.N,obj.N);
    for num = 1:size(C,2)
      bias_mat = bias_mat + x(1,num) * C{num};
      init_mat = init_mat + x0(1,num) * C{num};
    end
    Results{rep}.bias = diag(bias_mat);
    Results{rep}.init_bias = diag(init_mat);
    if isnan(T)
      Results{rep}.time = x(1,lC+1);
      Results{rep}.init_time = x0(1,lC+1);
    else
      Results{rep}.time = T;
      Results{rep}.init_time = T;
    end

    if v ~= 0
      % Setup trace plots
      if size(readout,2) == 0
        plot_time = [0:Results{rep}.time/1000:Results{rep}.time];
      else
        tt = Results{rep}.time + readout(1)/2;
        plot_time = [0:tt/1000:tt];
      end
      % natural evolution
      plot_nat = obj.trace(IN,OUT,plot_time);
      % evolution with control
      H = obj.H;
      for k = 2:length(ctrl.H)
        H = H + x(1,k-1) * C{k-1};
      end
      [V,e] = eig (H);
      e = diag(e);
      E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (plot_time), 'UniformOutput', false);
      plot_ctrl =  cellfun (@(x) abs(x(OUT, IN))^2, E);
      % Plot current result
      plot_result(rep, 2,4,1,5,'Current solution');
    end

    % Best solution
    if best == 0
      best = rep;
      if v ~= 0
        plot_result(best, 2,4,2,6,'Best solution');
      end
    else
      if Results{rep}.err < Results{best}.err
        best = rep;
        if v ~= 0
          plot_result(best, 2,4,2,6,'Best solution');
        end
      end
    end

    % Fastest solution
    if Results{rep}.err < min_err
      if fastest == 0
        fastest = rep;
        if v ~= 0
          plot_result(fastest,2,4,3,7,'Fastest solution');
        end
      elseif abs(Results{rep}.time - Results{fastest}.time) < eps(0)
        if Results{rep}.err < Results{fastest}.err
          fastest = rep;
          if v ~= 0
            plot_result(fastest,2,4,3,7,'Fastest solution');
          end
        end
      elseif Results{rep}.time < Results{fastest}.time
        fastest = rep;
        if v ~= 0
          plot_result(fastest,2,4,3,7,'Fastest solution');
        end
      end
    end

    if v ~= 0
      % Error vs time
      figure(1);
      subplot(2,4,4);
      hold on;
      if isnan(T)
        plot(log10(Results{rep}.err),x(end),'*b');
      else
        plot(log10(Results{rep}.err),T,'*b');
      end
      hold off;
      xlabel('log10(Error)');
      ylabel('Time');
      title('log10(Error) vs. Time');
      axis tight;

      % Error Histogram
      figure(1);
      subplot(2,4,8);
      Err(rep) = log10(Results{rep}.err);
      histogram(Err(1:rep),'FaceColor',[0 0 1],'FaceAlpha',1);
      title(sprintf('log10(Error) Histogram, %d runs', rep));
      axis tight;

      drawnow();
      refresh();
    end

  end

  if v ~= 0
    % Final plot
    clf

    if isnan(T)
      % Best solution
      if best > 0
        % Setup trace plots
        if size(readout,2) > 0
          tt = Results{best}.time + readout(1)/2;
        else
          tt = Results{best}.time;
        end
        plot_time = [0:tt/1000:tt];
        % natural evolution
        plot_nat = obj.trace(IN,OUT,plot_time);
        % evolution with control
        H = obj.H + diag(Results{best}.bias);
        [V,e] = eig (H);
        e = diag(e);
          E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (plot_time), 'UniformOutput', false);
        plot_ctrl =  cellfun (@(x) abs(x(OUT, IN))^2, E);
        plot_result(best, 3,3,1,4,'Best solution');
        subplot(3,2,5);
        plot_eigenstructure([V;e'],'best solution');
      end
      % Fastest solution
      if fastest > 0
        % Setup trace plots
        if size(readout,2) > 0
          tt = Results{fastest}.time + readout(1)/2;
        else
          tt = Results{fastest}.time;
        end
        plot_time = [0:tt/1000:tt];
        % natural evolution
        plot_nat = obj.trace(IN,OUT,plot_time);
        % evolution with control
        H = obj.H + diag(Results{fastest}.bias);
        [V,e] = eig (H);
        e = diag(e);
        E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (plot_time), 'UniformOutput', false);
        plot_ctrl =  cellfun (@(x) abs(x(OUT, IN))^2, E);
        plot_result(fastest, 3,3,2,5,'Fastest solution');
        subplot(3,2,6);
        plot_eigenstructure([V;e'],'fastest solution');
      end
      % Error vs time
      figure(1);
      subplot(3,3,3);
      times = arrayfun(@(x) (Results{x}.time),[1:size(Results,2)]);
      plot(log10(arrayfun(@(x) (Results{x}.err),[1:size(Results,2)])), times, '*b');
      xlabel('log(Error)');
      ylabel('Time');
      if isnan(T)
        ylim([max(0,min(times)) min(initT(1)*1.1,max(times))]);
      end
      title('log(Error) vs. Time');
      axis tight;
      % Error Histogram
      figure(1);
      subplot(3,3,6);
      Err = log10(arrayfun(@(x) (Results{x}.err),[1:size(Results,2)]));
      histogram(Err,'FaceColor',[0 0 1],'FaceAlpha',1);
      title(sprintf('log(Error) Histogram over %d runs', size(Results,2)));
      axis tight;
    else
      % Best solution
      if best > 0
        % Setup trace plots
        if size(readout,2) > 0
          tt = Results{best}.time + readout(1)/2;
        else
          tt = Results{best}.time;
        end
        plot_time = [0:tt/1000:tt];
        % natural evolution
        plot_nat = obj.trace(IN,OUT,plot_time);
        % evolution with control
        H = obj.H + diag(Results{best}.bias);
        [V,e] = eig (H);
        e = diag(e);
        E = cellfun (@(x) V * diag(exp(-i * x * e)) * V', num2cell (plot_time), 'UniformOutput', false);
        plot_ctrl =  cellfun (@(x) abs(x(OUT, IN))^2, E);
        plot_result(best, 3,2,1,3,'Best solution');
        subplot(3,2,[5 6]);
        plot_eigenstructure([V;e'],'best solution');
      end
      % Error Histogram
      figure(1);
      subplot(3,2,[2 4]);
      Err = log10(arrayfun(@(x) (Results{x}.err),[1:size(Results,2)]));
      histogram(Err,'FaceColor',[0 0 1],'FaceAlpha',1);
      title(sprintf('log(Error) Histogram over %d runs', size(Results,2)));
      axis tight;
    end

    drawnow();
    refresh();

    if size(readout,2) == 3
      figure(3);
      ctrl.test_static_st_dt_avg (x);
    end

    display(sprintf('Failed runs: %d\n', size(FailedRuns,2)));
  end

  % Plot results
  function plot_result (run, X, Y, trace_fig, bias_fig, str)
    % Plot traces
    figure(1);
    subplot(X,Y,trace_fig);
    plot(plot_time, plot_ctrl, plot_time, plot_nat);
    axis([0 plot_time(end) 0 1]);
    hold on;
    plot(plot_time(end),plot_ctrl(end),'*r');
    plot(plot_time(end),plot_nat(end),'*g');
    hold off;
    %legend('control','no control', 'Location', 'best')
    xlabel('time in 1/J')
    ylabel(sprintf('transition probability p_{%i%i}',IN,OUT))
    title(sprintf('%s, N=%d, e=%.6g', str, obj.N, Results{run}.err));
    axis tight;
    % Plot bias
    figure(1);
    subplot(X,Y,bias_fig);
    bar(1:obj.N,Results{run}.bias,'b');
    hold on;
    bar(1:obj.N,Results{run}.init_bias,0.2,'g');
    hold off;
    xlabel('Spin #');
    ylabel('Bias');
    title(sprintf('Bias |%d>-|%d>, T=%.6g, T0=%.6g',IN,OUT,Results{run}.time,Results{run}.init_time));
    axis tight;
  end

  function plot_eigenstructure(V,str)
    colormap(flipud(gray));
    VV = abs(V);
    clim = [min(min(VV(1:end-1,1:end))) max(max(VV(1:end-1,1:end)))];
    im = imagesc(VV,clim);
    imAxes =get(im,'parent');
    hFig = get(imAxes,'parent');
    fs = getBestFontSize(imAxes);
    axis off;
    [rows, cols] = size(V);
    midValue = mean(get(gca,'CLim'));
    ci = (VV < midValue) + 1;
    cmap = colormap ();
    [mx my] = size(cmap);
    cmap = [cmap(1,:); cmap(mx,:)];
    textHandles = zeros(size(V))';
    for i = 1:rows
      for j = 1:cols
        c = cmap(ci(i,j),:);
        if i == IN
          c = [0 ((2*c(1)-1)/4+3/4) 0];
        elseif i == OUT
          c = [((2*c(1)-1)/4+3/4) 0 0];
        elseif i == rows
          c = [1 1 1];
        end
        textHandles(j,i) = text(j,i,num2str(V(i,j),'%.3g'),...
                                'color', c,...
                                'horizontalAlignment','center','verticalAlignment','middle',...
                                'fontsize',fs,'clipping','on','visible','on');
      end
    end
    set(imAxes,'UserData',textHandles);
    title(['Eigenstructure of ' str]);
  end

  function fs = getBestFontSize(imAxes)
    % Try to keep font size reasonable for text
    hFig = get(imAxes,'Parent');
    magicNumber = 80;
    nrows = diff(get(imAxes,'YLim'));
    ncols = diff(get(imAxes,'XLim'));
    if ncols < magicNumber && nrows < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/max(nrows,ncols);
    elseif ncols < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/ncols;
    elseif nrows < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/nrows;
    else
      ratio = 1;
    end
    fs = min(10,ceil(ratio/5));    % the gold formula
  end

end
