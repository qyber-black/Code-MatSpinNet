function [F,dF] = dexpma(M,dM,v)
  % [F,dF] = dexpma(M,dM,v)
  %
  % Gradient helper for eal_static_st

  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Frank C Langbein <frank@langbein.org>, Cardiff University
  % SPDX-FileCopyrightText: Copyright (C) 2011-2019, 2022 Sophie M Shermer <lw1660@gmail.com>, Swansea University
  % SPDX-License-Identifier: AGPL-3.0-or-later
  
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