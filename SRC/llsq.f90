subroutine llsq0 ( n, x, y, a )

!*****************************************************************************80
!
!! LLSQ solves a linear least squares problem matching y=a*x to data.
!
!  Discussion:
!
!    A formula for a line of the form Y = A * X is sought, which
!    will minimize the root-mean-square error to N data points ( X(I), Y(I) );
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2019
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the data points.
!
!    Output, real ( kind = 8 ) A the slope of the 
!    least-squares approximant to the data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) bot
  real ( kind = 8 ) top
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) y(n)
!
!  Special case.
!
  if ( n == 1 ) then
    if ( x(1) == 0.0D+00 ) then
      a = 1.0D+00
    else
      a = y(1) / x(1)
    end if
    return
  end if

  top = dot_product ( x(1:n), y(1:n) )
  bot = dot_product ( x(1:n), x(1:n) )

  a = top / bot

  return
end
subroutine llsq ( n, x, y, a, b )

!*****************************************************************************80
!
!! LLSQ solves a linear least squares problem matching a line to data.
!
!  Discussion:
!
!    A formula for a line of the form Y = A * X + B is sought, which
!    will minimize the root-mean-square error to N data points ( X(I), Y(I) );
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) X(N), Y(N), the coordinates of the data points.
!
!    Output, real ( kind = 8 ) A, B, the slope and Y-intercept of the 
!    least-squares approximant to the data.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bot
  real ( kind = 8 ) top
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xbar
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ybar
!
!  Special case.
!
  if ( n == 1 ) then
    a = 0.0D+00
    b = y(1)
    return
  end if
!
!  Average X and Y.
!
  xbar = sum ( x(1:n) ) / real ( n, kind = 8 )
  ybar = sum ( y(1:n) ) / real ( n, kind = 8 )
!
!  Compute Beta.
!
  top = dot_product ( x(1:n) - xbar, y(1:n) - ybar )
  bot = dot_product ( x(1:n) - xbar, x(1:n) - xbar )

  a = top / bot

  b = ybar - a * xbar

  return
end