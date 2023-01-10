function step = hammersley_step_get ( )

%*****************************************************************************80
%
%% HAMMERSLEY_STEP_GET gets the "step" of the leaped Hammersley subsequence.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 July 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Output, integer STEP, the step of the leaped Hammersley subsequence.
%

% SPDX-FileCopyrightText: John Burkardt, 2004
% SPDX-License-Identifier: LGPL-3.0-only

  global hammersley_STEP

  step = hammersley_STEP;

  return
end
