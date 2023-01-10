function example_hes
  % QSN package example
  %
  % Higher excitation subspace example

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct hes chain and display
  disp('=======================================================================');
  disp('Constructing an XX chain of length 5');
  N = 5;
  chain = qsn.QSN('chain', N)
  key ();

  % Hes chain
  disp('Constructing higher excitation subspace 0:3 for chain');
  J = chain.H;
  Z = zeros(N,N);
  M = ceil(N/2);
  hes = qsn.QSN('hes',{J J Z},0:M)
  key ();

  % Plot hes chain
  disp('Plotting higher excitation subspace 0:3 for chain')
  figure(1);
  clf;
  hes.plot ();
  key ();

  function key ()
    disp('Press key to continue');
    pause;
  end

end