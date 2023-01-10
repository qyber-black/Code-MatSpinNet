function [N,x] = opt_attainable (obj,ispin,jspin,x,s)
  % N = obj . opt_attainable (ispin, jspin, x, s, shortest) - Target function for finding attainability
  %
  %   ispin - Initial spin
  %   jspin - Target spin
  %   x     - non-uniform scaling vector
  %   s     - scalar scaling parameter
  %   obj   - Quantum spin network object
  %   N     - Number of fulfilled even/odd conditions

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later

  [~,q,~,~,~,info] = obj.is_attainable(ispin,jspin, [x s]);
  if q == 0
    N = 100;
  else
    N = size(info,1) - sum(info(1:end,3));
  end
end