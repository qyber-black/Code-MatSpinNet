function hammersley_step_set ( step )

%*****************************************************************************80
%
%% HAMMERSLEY_STEP_SET sets the step of the leaped Hammersley subsequence.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    29 January 2005
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer STEP, the step of the leaped Hammersley subsequence.
%

% SPDX-FileCopyrightText: John Burkardt, 2004
% SPDX-License-Identifier: LGPL-3.0-only

  global hammersley_STEP

  step = floor ( step );

  if ( ~halham_step_check ( step ) )
    error ( 'HAMMERSLEY_STEP_SET - Fatal error!' );
  end

  hammersley_STEP = step;

  return
end
