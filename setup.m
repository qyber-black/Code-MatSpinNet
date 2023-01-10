function setup(a)
  % Setup QSN code
  %  setup('build') - compile code
  %  setup()        - set paths (or any other argument)

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 SM Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  if exist('a','var') && strcmp(a,'build')
    % Compile KDE code
    cd '@kde/mex'
    makemex
    cd '../..'
  end
  addpath '.' './GA' 
end
