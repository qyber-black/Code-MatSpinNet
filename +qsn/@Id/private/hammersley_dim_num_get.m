function dim_num = hammersley_dim_num_get ( )

%*****************************************************************************80
%
%% HAMMERSLEY_DIM_NUM_GET gets the dimension of the leaped Hammersley subsequence.
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
%    Output, integer DIM_NUM, the current value of the dimension.   
%

% SPDX-FileCopyrightText: John Burkardt, 2004
% SPDX-License-Identifier: LGPL-3.0-only

  global hammersley_DIM_NUM

  dim_num = hammersley_DIM_NUM;

  return
end
