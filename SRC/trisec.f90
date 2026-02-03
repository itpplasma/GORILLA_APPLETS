module trisec_mod
!
! Module for triangle-line intersection calculations.
!
! Contains routines for:
!   - Finding intersections of a line with a triangulated surface
!   - Computing line-plane intersections
!   - Testing if a point lies inside a triangle
!
  use, intrinsic :: iso_fortran_env, only: dp => real64

  implicit none

  private
  public :: find_surface_intersection, line_plane_intersection, point_in_triangle
  public :: trisec  ! Legacy name, deprecated

contains

! ====================================================================
subroutine find_surface_intersection(start_point, direction, n_triangles, &
                                     vertices1, vertices2, vertices3, &
                                     intersection_point, min_distance, hit_triangle)
!
! Finds the nearest intersection of a ray with a triangulated surface.
!
! A ray is cast from start_point in the given direction. The routine
! searches through all triangles and returns the nearest intersection
! point (in the positive direction).
!
! Input:
!   start_point(3)    - Origin point of the ray
!   direction(3)      - Direction vector of the ray
!   n_triangles       - Number of triangles in the surface
!   vertices1(3,n)    - First vertex of each triangle
!   vertices2(3,n)    - Second vertex of each triangle
!   vertices3(3,n)    - Third vertex of each triangle
!
! Output:
!   intersection_point(3) - Coordinates of the nearest intersection
!   min_distance          - Distance from start_point to intersection
!   hit_triangle          - Index of the intersected triangle (0 if none)
!
  real(dp), intent(in)  :: start_point(3)
  real(dp), intent(in)  :: direction(3)
  integer,  intent(in)  :: n_triangles
  real(dp), intent(in)  :: vertices1(3,n_triangles)
  real(dp), intent(in)  :: vertices2(3,n_triangles)
  real(dp), intent(in)  :: vertices3(3,n_triangles)
  real(dp), intent(out) :: intersection_point(3)
  real(dp), intent(out) :: min_distance
  integer,  intent(out) :: hit_triangle

  real(dp) :: edge21(3), edge31(3), normal(3), plane_offset
  real(dp) :: candidate_point(3), distance
  integer  :: i, inside_flag

  hit_triangle = 0
  min_distance = huge(min_distance)

  do i = 1, n_triangles
     ! Compute triangle edges
     edge21 = vertices2(:,i) - vertices1(:,i)
     edge31 = vertices3(:,i) - vertices1(:,i)

     ! Compute plane normal via cross product: edge21 x edge31
     normal(1) = edge21(2)*edge31(3) - edge31(2)*edge21(3)
     normal(2) = edge21(3)*edge31(1) - edge31(3)*edge21(1)
     normal(3) = edge21(1)*edge31(2) - edge31(1)*edge21(2)

     ! Plane offset: D = -normal dot vertex1
     plane_offset = -dot_product(normal, vertices1(:,i))

     ! Find intersection with the triangle's plane
     call line_plane_intersection(normal, plane_offset, direction, start_point, distance)

     ! Check if intersection is in positive direction and closer than current best
     if (distance > 0.0d0 .and. distance < min_distance) then
        ! Compute candidate intersection point
        candidate_point = start_point + distance * direction

        ! Test if point lies inside the triangle
        inside_flag = sign(1.0d0, point_in_triangle(vertices1(:,i), &
                           vertices2(:,i), vertices3(:,i), candidate_point))

        if (inside_flag > 0) then
           hit_triangle = i
           min_distance = distance
           intersection_point = candidate_point
        endif
     endif
  enddo

end subroutine find_surface_intersection

! ====================================================================
subroutine trisec(x0,y0,z0,delx,dely,delz,n,x1,y1,z1,x2,y2,z2,x3,y3,z3,  &
                  xps,yps,zps,rhomn,ntri)
!
! Legacy wrapper for find_surface_intersection.
! DEPRECATED: Use find_surface_intersection instead.
!
! Calculates the intersection of a straight line starting at (x0,y0,z0)
! in direction (delx,dely,delz) with a surface of N triangles.
!
  real(dp), intent(in)  :: x0, y0, z0
  real(dp), intent(in)  :: delx, dely, delz
  integer,  intent(in)  :: n
  real(dp), intent(in)  :: x1(n), y1(n), z1(n)
  real(dp), intent(in)  :: x2(n), y2(n), z2(n)
  real(dp), intent(in)  :: x3(n), y3(n), z3(n)
  real(dp), intent(out) :: xps, yps, zps
  real(dp), intent(out) :: rhomn
  integer,  intent(out) :: ntri

  real(dp) :: start_point(3), direction(3), intersection_point(3)
  real(dp) :: vertices1(3,n), vertices2(3,n), vertices3(3,n)
  integer  :: i

  ! Pack scalar arguments into vectors
  start_point = [x0, y0, z0]
  direction = [delx, dely, delz]

  ! Repack vertex arrays from (n) to (3,n) layout
  do i = 1, n
     vertices1(:,i) = [x1(i), y1(i), z1(i)]
     vertices2(:,i) = [x2(i), y2(i), z2(i)]
     vertices3(:,i) = [x3(i), y3(i), z3(i)]
  enddo

  ! Call the new implementation
  call find_surface_intersection(start_point, direction, n, &
                                 vertices1, vertices2, vertices3, &
                                 intersection_point, rhomn, ntri)

  ! Unpack results
  xps = intersection_point(1)
  yps = intersection_point(2)
  zps = intersection_point(3)

end subroutine trisec

! ====================================================================
subroutine line_plane_intersection(normal_vec, plane_offset, direction, point, rho)
!
! Calculates the intersection of a straight line with a plane.
!
! The plane is defined by: normal_vec dot r + plane_offset = 0
!   where normal_vec is the plane normal vector
!   and plane_offset = -(normal_vec dot r_plane) for any point r_plane on the plane
!
! The line is defined by: r(t) = point + t * direction
!   where point is a point on the line
!   and direction is the unit direction vector
!
! Output: rho = distance along the line from 'point' to the intersection
!         (positive if intersection is in the direction of travel)
!
! If the line is parallel to the plane (no intersection), rho = huge(rho)
!
  real(dp), intent(in)  :: normal_vec(3)   ! Plane normal vector
  real(dp), intent(in)  :: plane_offset    ! Plane offset D in equation: normal dot r + D = 0
  real(dp), intent(in)  :: direction(3)    ! Line direction (unit vector)
  real(dp), intent(in)  :: point(3)        ! Point on the line
  real(dp), intent(out) :: rho             ! Distance to intersection

  real(dp) :: normal_dot_direction, normal_dot_point_plus_offset
  real(dp), parameter :: eps = 1.0d-7

  ! Dot product of normal with line direction
  ! If zero, line is parallel to plane (no intersection)
  normal_dot_direction = dot_product(normal_vec, direction)

  ! Evaluate plane equation at the point: normal dot point + plane_offset
  ! This measures signed distance of point from plane (scaled by |normal|)
  normal_dot_point_plus_offset = dot_product(normal_vec, point) + plane_offset

  ! Check if line is parallel to plane
  if (abs(normal_dot_direction) <= abs(normal_dot_point_plus_offset) * eps) then
     rho = huge(rho)  ! No intersection (line parallel to plane)
  else
     ! Solve for rho: at intersection, normal dot (point + rho*direction) + plane_offset = 0
     ! => rho = -(normal dot point + plane_offset) / (normal dot direction)
     rho = -normal_dot_point_plus_offset / normal_dot_direction
  endif

end subroutine line_plane_intersection

! ====================================================================
function point_in_triangle(vertex1, vertex2, vertex3, test_point)
!
! Determines if a point lies inside a triangle (assuming coplanar).
!
! The test uses the sign of cross products: if the point is inside the triangle,
! all cross products (point-to-vertex1) x (point-to-vertex2),
! (point-to-vertex2) x (point-to-vertex3), and (point-to-vertex3) x (point-to-vertex1)
! will point in the same direction (same sign when dotted together).
!
! Input:
!   vertex1, vertex2, vertex3 - the three vertices of the triangle
!   test_point - the point to test
!
! Output:
!   point_in_triangle > 0 if test_point is inside the triangle
!   point_in_triangle < 0 if test_point is outside the triangle
!
  real(dp), intent(in) :: vertex1(3), vertex2(3), vertex3(3)
  real(dp), intent(in) :: test_point(3)
  real(dp) :: point_in_triangle

  real(dp) :: vec_to_v1(3), vec_to_v2(3), vec_to_v3(3)
  real(dp) :: cross_12(3), cross_23(3), cross_31(3)
  real(dp) :: dot_12_23, dot_31_23

  ! Build vectors from test point to each vertex
  vec_to_v1 = vertex1 - test_point
  vec_to_v2 = vertex2 - test_point
  vec_to_v3 = vertex3 - test_point

  ! Compute cross products of consecutive vectors
  ! cross_12 = vec_to_v1 x vec_to_v2
  cross_12(1) = vec_to_v1(2)*vec_to_v2(3) - vec_to_v1(3)*vec_to_v2(2)
  cross_12(2) = vec_to_v1(3)*vec_to_v2(1) - vec_to_v1(1)*vec_to_v2(3)
  cross_12(3) = vec_to_v1(1)*vec_to_v2(2) - vec_to_v1(2)*vec_to_v2(1)

  ! cross_23 = vec_to_v2 x vec_to_v3
  cross_23(1) = vec_to_v2(2)*vec_to_v3(3) - vec_to_v3(2)*vec_to_v2(3)
  cross_23(2) = vec_to_v2(3)*vec_to_v3(1) - vec_to_v3(3)*vec_to_v2(1)
  cross_23(3) = vec_to_v2(1)*vec_to_v3(2) - vec_to_v3(1)*vec_to_v2(2)

  ! Check if cross_12 and cross_23 point in the same direction
  dot_12_23 = dot_product(cross_12, cross_23)
  point_in_triangle = dot_12_23
  if (dot_12_23 < 0.0d0) return  ! Already outside, no need to check further

  ! cross_31 = vec_to_v3 x vec_to_v1
  cross_31(1) = vec_to_v3(2)*vec_to_v1(3) - vec_to_v3(3)*vec_to_v1(2)
  cross_31(2) = vec_to_v3(3)*vec_to_v1(1) - vec_to_v3(1)*vec_to_v1(3)
  cross_31(3) = vec_to_v3(1)*vec_to_v1(2) - vec_to_v3(2)*vec_to_v1(1)

  ! Check if cross_31 and cross_23 point in the same direction
  dot_31_23 = dot_product(cross_31, cross_23)

  ! Point is inside only if all cross products point in the same direction
  point_in_triangle = min(dot_12_23, dot_31_23)

end function point_in_triangle

end module trisec_mod
