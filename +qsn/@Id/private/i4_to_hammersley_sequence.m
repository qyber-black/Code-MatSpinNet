function r = i4_to_hammersley_sequence ( dim_num, n, step, seed, leap, base )

%*****************************************************************************80
%
%% I4_TO_HAMMERSLEY_SEQUENCE: next N elements of an DIM_NUM-dimensional Hammersley sequence.
%
%  Discussion:
%
%    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
%    sequences, each generated by a particular base.  If the base is 
%    greater than 1, a standard 1-dimensional
%    van der Corput sequence is generated.  But if the base is 
%    negative, this is a signal that the much simpler sequence J/(-BASE) 
%    is to be generated.  For the standard Hammersley sequence, the
%    first spatial coordinate uses a base of (-N), and subsequent
%    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
%    This program allows the user to specify any combination of bases,
%    included nonprimes and repeated values.
%
%    This routine selects elements of a "leaped" subsequence of the
%    Hammersley sequence.  The subsequence elements are indexed by a
%    quantity called STEP, which starts at 0.  The STEP-th subsequence
%    element is simply element
%
%      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
%
%    of the original Hammersley sequence.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    05 May 2008
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    John Hammersley,
%    Monte Carlo methods for solving multivariable problems,
%    Proceedings of the New York Academy of Science,
%    Volume 86, 1960, pages 844-874.
%
%    Ladislav Kocis, William Whiten,
%    Computational Investigations of Low-Discrepancy Sequences,
%    ACM Transactions on Mathematical Software,
%    Volume 23, Number 2, 1997, pages 266-294.
%
%  Parameters:
%
%    Input, integer DIM_NUM, the spatial dimension.
%    1 <= DIM_NUM is required.
%
%    Input, integer N, the number of elements of the sequence.
%
%    Input, integer STEP, the index of the subsequence element.
%    0 <= STEP is required.
%
%    Input, integer SEED(DIM_NUM), the sequence index corresponding
%    to STEP = 0.
%
%    Input, integer LEAP(DIM_NUM), the succesive jumps in the sequence.
%
%    Input, integer BASE(DIM_NUM), the bases.
%
%    Output, real R(DIM_NUM,N), the next N elements of the
%    leaped subsequence, beginning with element STEP.
%

% SPDX-FileCopyrightText: John Burkardt, 2008
% SPDX-License-Identifier: LGPL-3.0-only

  fiddle = 0;
  dim_num = floor ( dim_num );
  n = floor ( n );
  step = floor ( step );
  seed(1:dim_num) = floor ( seed(1:dim_num) );
  leap(1:dim_num) = floor ( leap(1:dim_num) );
  base(1:dim_num) = floor ( base(1:dim_num) );
%
%  Check the input.
%
  if ( ~halham_dim_num_check ( dim_num ) )
    error ( 'I4_TO_HAMMERSLEY_SEQUENCE - Fatal error!' );
  end

  if ( ~halham_n_check ( n ) )
    error ( 'I4_TO_HAMMERSLEY_SEQUENCE - Fatal error!' );
  end

  if ( ~halham_step_check ( step ) )
    error ( 'I4_TO_HAMMERSLEY_SEQUENCE - Fatal error!' );
  end

  if ( ~halham_seed_check ( dim_num, seed ) )
    error ( 'I4_TO_HAMMERSLEY_SEQUENCE - Fatal error!' );
  end

  if ( ~halham_leap_check ( dim_num, leap ) )
    error ( 'I4_TO_HAMMERSLEY_SEQUENCE - Fatal error!' );
  end

  if ( ~hammersley_base_check ( dim_num, base ) )
    error ( 'I4_TO_HAMMERSLEY_SEQUENCE - Fatal error!' );
  end
%
%  Calculate the data.
%
  r(1:dim_num,1:n) = 0.0;
  
  for i = 1: dim_num

    if ( 1 < base(i) )

      seed2(1:n) = seed(i) + step * leap(i) : leap(i) : ...
                   seed(i) + ( step + n - 1 ) * leap(i);

      base_inv = 1.0 / base(i);
  
      while ( any ( seed2 ~= 0 ) )
        digit(1:n) = mod ( seed2(1:n), base(i) );
        r(i,1:n) = r(i,1:n) + digit(1:n) * base_inv;
        base_inv = base_inv / base(i);
        seed2(1:n) = floor ( seed2(1:n) / base(i) );
      end
%
%  In the following computation, the value of FIDDLE can be:
%
%    0,   for the sequence 0/N, 1/N, ..., N-1/N
%    1,   for the sequence 1/N, 2/N, ..., N/N
%    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
%
    else

      temp(1:n) = seed(i) + step * leap(i) : leap(i) : ...
                 seed(i) + ( step + n - 1 ) * leap(i);

      r(i,1:n) = ( mod ( temp(1:n), ( -base(i) + 1 ) ) + fiddle ) / ( -base(i) );

    end

  end

  return
end
