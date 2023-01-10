function example_networks
  % QSN package example
  %
  % Analyse geometry of a quantum spin ring range to demonstrate basic
  % usage of qsn.QSNRange and qsn.DSRange classes.

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct rings and display
  disp('=======================================================================');
  disp('Construct range of spin network rings with nodes 11 to 31');
  nets = qsn.QSNRange (@(n) qsn.QSN('ring', n), 'Nodes', 11:31)
  key ();

  disp('Construct distance space range for this network range');
  dists = qsn.DSRange (nets)
  key ();

  disp ('Plot maximum Gromov deltas for distance space range');
  [delta4,scaled_delta4] = dists.max_Gromov4pt (1) 
  key ();

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct engineered rings and display
  disp('=======================================================================');
  disp('Construct range of spin network rings with nodes 11 to 31 and');
  disp('potential in the middle of 0 100 500 1000');
  nets = qsn.QSNRange (@(n,z) qsn.QSN('ring', n, 'XX', [zeros(1,floor(n/2)) z]), 'Nodes', 11:31, 'Zeta', [0 100 500 1000])
  key ();

  disp('Construct distance space range for this network range');
  dists = qsn.DSRange (nets)
  key ();

  disp ('Plot maximum Gromov deltas for distance space range in 2D');
  [delta4,scaled_delta4] = dists.max_Gromov4pt (1) 
  key ();

  disp ('Plot maximum Gromov deltas for distance space range as curves');
  [delta4,scaled_delta4] = dists.max_Gromov4pt (1, true, true) 
  key ();

  function key () 
    disp('Press key to continue');
    pause;
  end

end