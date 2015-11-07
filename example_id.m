function example_control
  % QSN package example
  %
  % Show basic usage of qsn.Id.
  %

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct a ring network
  disp('=======================================================================');
  disp('Construct an XX ring network with 6 nodes');
  net = qsn.QSN ('ring', ones(1,8)*1.2)
  key ();

  % Setup identification problem
  disp('Setup identification problem');
  id = qsn.Id (net, [5 15; .5 1.5], 50)
  key ();

  % Run identification
  disp('Finding parameters...');
  [H,err,D] = id.find (1, [0 40], 10, 10, 10);
  key ();

  function key ()
    disp('Press key to continue');
    pause;
  end

end
