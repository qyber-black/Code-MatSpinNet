function example_control
  % QSN package example
  %
  % Show basic usage of qsn.Control and qsn.Opt by controlling an XX chain with 11 nodes.
  
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct a chain network
  disp('=======================================================================');
  disp('Construct an XX chain network with 11 nodes');
  chain = qsn.QSN ('chain', 11)
  key ();

  % Setup switch control problem
  disp('Construct an XX spin network chain with 11 nodes');
  C = zeros (11, 11);
  C(1,1) = 1000;
  ctrl = qsn.Control (chain, 'state_transfer', C,  1, 11, 10, NaN)
  key ();

  % Solve control problem
  disp('Solving...');
  [x,err,exec_time,exit_flag,output,x0] = ctrl.solve_switch (qsn.Opt ('fmin'), 10)
  key ();

  function key ()
    disp('Press key to continue');
    pause;
  end

end