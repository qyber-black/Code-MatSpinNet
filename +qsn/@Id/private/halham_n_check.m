function value = halham_n_check ( n )

%*****************************************************************************80
%
%% HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 September 2004
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer N, the number of points in the subsequence.
%
%    Output, logical VALUE, is true if N is legal.
%

% SPDX-FileCopyrightText: John Burkardt, 2004
% SPDX-License-Identifier: LGPL-3.0-only

  if ( n < 1 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'HALHAM_N_CHECK - Fatal error!\n' );
    fprintf ( 1, '  Input value of N < 1!\n' );
    fprintf ( 1, '  N = %d\n', n );
    value = 0;
  else
    value = 1;
  end

  return
end
