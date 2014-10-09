# GeometricalPredicates

[![Build Status](https://travis-ci.org/skariel/GeometricalPredicates.jl.svg?branch=master)](https://travis-ci.org/skariel/GeometricalPredicates.jl)

Fast, robust 2D and 3D geometrical predicates on generic point types.
Implementation follows algorithms described in the [Arepo paper](http://arxiv.org/abs/0901.4107)
and used (for e.g.) in the [Illustris Simulation](http://www.illustris-project.org/). License: MIT. Bug reports welcome!

How does it work?
--------------------
Calculations initially performed on Float64 while bounding max
absolute errors. If unable to determine result fall back to exact
calculation using BigInts. This is a form of floating point filtering.
Most calculations are cached for fast repeated testing of
incircle/intriangle predicates

How to use?
--------------
###Installation
```Julia
Pkg.clone("git@github.com:skariel/GeometricalPredicates.jl.git")
```
###Points
```Julia
using GeometricalPredicates

# create a 2D point in (x, y) = (1.1, 1.9)
mypoint = Point(1.1, 1.9)
typeof(mypoint) # -> Point2D

# create a 3D point in (x, y, z) = (1.1, 1.9, 1.5)
mypoint = Point(1.1, 1.9, 1.5)
typeof(mypoint) # -> Point3D

# getting coordinates:
getx(mypoint) # -> 1.1
gety(mypoint) # -> 1.9
getz(mypoint) # -> 1.5
```
`Point2D` inherits from `AbstractPoint2D`and `Point3D` inherits from `AbstractPoint3D`.
You can implement custom point types by inheriting from these abstract types. These
custom point types can be used with the rest of the package:
```
type MyCustomPointType <: AbstractPoint2D
    _x::FLoat64
    _y::Float64
    _mass::Float64
end

getx(p::MyCustomPointType) = p._x
gety(p::MyCustomPointType) = p._y
```
implementing `getx`, `gety`, and `getz` for 3D points is necessary
as this is the interface the package is expecting.

The point coordinates must reside in a region `1.0 <= x < 2.0`. Read below on
why is this limitation necessary. For convenience there are 2 constants defined,
`min_coord` and `max coord` representing the minimal and maximal feasible values
of coordinates.

###Triangles and Tetrahedrons (..aka Primitives)



Usage notes
-------------
- all point types must implement getx(p), gety(p) and (if 3D point) getz(p).
  These methods should return Float64. The predicated are generic in the sense
  that they accept any such point carrying additional user defined data.
- all point coordinates must be in the range [1, 2), i.e. including 1.0 and
  ending with inclusion of at 2.0-e where thiss represents the Float64 number
  which is closest to 2.0
- beware of creating points (2D or 3D) using somthing like Point2D(3, 4) -
  it will return the point belonging to a Peano-Hilbert key `3` in
  an 2^4 X 2^4 grid. If you want to create a point located at 3, 4 use
  Point2D(3.0, 4.0)

`incircle` method: determines if the point 
      lies inside (ret->1), exactly on (ret->0),
      outside (ret->-1), or NA (ret->2) of a circle or sphere
      defined by the primitive

`intriangle` method: determines if the point
      lies inside (ret->1) or outside (ret->negative num) the primitive.
      The negative number is the negative index ot the triangle point
      which test point is in front of. If test point is in front of two
      triangle points then one is chosen in an undefined mannar.
      If point is exactly on one of the sides it returns the
      index+1 of the point that is in front of said side
      i.e. for a point infront of Triangle.a return 2,
           for a point infront of Triangle.b return 3, etc.
      If test point is exactly on one of the triangle points then one
      side is chosen in an undefined mannar.

`peanokey` method: this is the scale-dependent Peano-Hilbert interface.
      returns the Peano-Hilbert key for the given
      point in an nXn grid starting with 0 at (1.0, 2.0-e) and
      ending at nXn-1 at (2.0-e,2.0-e) for 2D, and starting at
      (1.0, 1.0, 1.0) with zero and ending with nXnXn-1 at
      (1.0, 1.0, 2.0-e) for 3D. `n` here is 2^bits, `bits` being a parameter
      passed to the function. For inverse Peano-Hilbert keys just create
      points using integers such as Point3D(1, 5) this will create a 3D point
      with peano index of 3 on a grid of size 2^5 X 2^5 X 2^5

