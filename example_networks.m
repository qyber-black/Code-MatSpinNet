function example_networks
  % QSN package example
  %
  % Analyse geometry of a quantum spin chains and rings to demonstrate
  % basic usage of qsn.QSN and qsn.DS classes.
  % 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct chain and display
  disp('=======================================================================');
  disp('Constructing an XX chain of length 5');
  chain = qsn.QSN('chain', 5)
  key ();

  % Chain figure
  disp('Plotting information about chain');
  figure(1);
  clf;
  chain.plot ();
  key ();

  % Estimate maximal probability within time interval
  disp ('Estimating probabilities attained between 0 and 10');
  [prob_max,time] = chain.prob_max (0, 10, 5, 6, 1)
  key ();

  % Plot dynamics of chain
  clf;
  disp ('Plotting dynamics of chain from 0 to 10');
  chain.plot_dynamics(.5, 0:.05:10);
  key ();

  disp ('Plotting ddynamics of chain from 0 to 10 for node 1');
  chain.plot_dynamics(.5, 0:.05:10, -1);
  key ();

  % Construct distance space
  disp ('Construct distances from probabilities');
  dist = qsn.DS(chain)
  dist.plot ();
  key ();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct engineered chain and display
  disp('=======================================================================');
  disp('Constructing an engineered Heisenberg chain of length 5');
  chain = qsn.QSN('chain', 5, 'H', [ 0 0 1000 0 0 ])
  key ();

  % Chain figure
  disp('Plotting information about chain');
  figure(1);
  clf;
  chain.plot ();
  key ();

  % Construct distance space
  disp ('Construct distances from probabilities');
  dist = qsn.DS(chain)
  dist.plot ();
  key ();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct ring and display
  disp('=======================================================================');
  disp('Constructing an XX ring of length 5');
  ring = qsn.QSN('ring', 5)
  key ();

  % Chain figure
  disp('Plotting information about ring');
  figure(1);
  clf;
  ring.plot ();
  key ();

  % Estimate maximal probability within time interval
  disp ('Estimating probabilities attained between 0 and 10');
  [prob_max,time] = ring.prob_max (0, 10, 5, 6, 1)
  key ();

  % Plot dynamics of ring
  clf;
  disp ('Plotting dynamics of ring from 0 to 10');
  ring.plot_dynamics(.5, 0:.05:10);
  key ();

  disp ('Plotting ddynamics of ring from 0 to 10 for node 1');
  ring.plot_dynamics(.5, 0:.05:10, -1);
  key ();

  % Construct distance space
  disp ('Construct distances from probabilities');
  dist = qsn.DS(ring)
  dist.plot ();
  key ();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct engineered ring and display
  disp('=======================================================================');
  disp('Constructing an engineered Heisenberg ring of length 13');
  ring = qsn.QSN('ring', 13, 'H', [0 0 0 0 0 0 1000 0 0 0 0 0 0])
  key ();

  % Chain figure
  disp('Plotting information about ring');
  figure(1);
  clf;
  ring.plot ();
  key ();

  % Construct distance space
  disp ('Construct distances from probabilities');
  dist = qsn.DS(ring)
  dist.plot ();
  key ();

  disp ('Checking metric:');
  dist.check_metric (true);
  key ();

  disp ('Fixing metric:');
  M = dist.generate_metric ('mean', 'mean', 'isoceles', @(x) mean(x), true);
  key ();
  
  disp ('Displaying resulting metric:')
  disp(M);
  M.plot ();
  M.check_metric (true);
  key ();
  
  function key () 
    disp('Press key to continue');
    pause;
  end

end