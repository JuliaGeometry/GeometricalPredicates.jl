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
```Julia
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
A triangle is the 2D primitive, and a tetrahedron is the 3D primitive
```Julia
# create a triangle using 3 points
a = Point(1.1, 1.1)
b = Point(1.9, 1.1)
c = Point(1.1, 1.9)
mytriangle = Primitive(a, b, c)
typeof(mytriangle) # -> UnOrientedTriangle{Point2D}
```
The triangle is unoriented in the sense that orientation is not-known in advance,
it is not immutable and it could change if points in the triangle are updated.
The orientation needs to be calculated when the triangle is created and when
points within are updated. Read below for the definition of orientation.nThe
triangle could be created using any points inheriting from `AbstractPoint2D`
which impement `getx` and `gety`, or using coordinates directly:
```Julia
mytriangle = Primitive(1.1, 1.1, 1.9, 1.1, 1.1, 1.9)

# Getting point `a` in the triangle
geta(mytriangle) # -> Point2D(1.1, 1.1)
getb(mytriangle) # -> Point2D(1.9, 1.1)
getc(mytriangle) # -> Point2D(1.1, 1.9)
```
The same goes for tetrahedrons, except we now use 4 3D points instead of 3 2D ones:
```Julia
# create a tetrahedron using 4 points
a = Point(1.1, 1.1, 1.1)
b = Point(1.9, 1.1, 1.1)
c = Point(1.1, 1.9, 1.1)
d = Point(1.1, 1.1, 1.9)
mytetraedron = Primitive(a, b, c, d)\
typeof(mytetrahedron) # -> UnOrientedTetrahedron{Point3D}
```
For certain applications we use primitives with known orientation, in those cases
there should be no need to calculate it. This is achieved in this package
by passing an `orientation` flag to `Primitive` creation function:
```Julia
mytetrahedron = Primitive(a, b, c, d, positivelyoriented)
typeof(mytetrahedron) # -> PositivelyOrientedTetrahedron{Point3D}
orientation(mytetrahedron) # -> constant 1, not calculated
mytetrahedron = Primitive(a, b, c, d, negativelyoriented)
typeof(mytetrahedron) # -> NegativelyOrientedTetrahedron{Point3D}
orientation(mytetrahedron) # -> constant -1, not calculated
```
Note that whe the primitive is oriented the real orientation is never calculated.
It is assumed that the user knows what he's doing. If in doubt, just use unoriented
primitives at the cost of actual calculation.

Updtating points in primitives can be done with `seta`, `setb`, etc. methods:
```Julia
seta(mytriangle, Point(1.7, 1.7))
```
Updating a point in a primitive will fire all relevant pre-calculations. I.e. if the triangle
is unoriented then orientation will be calculated. If is oriented then still a few other
pre-calculations will be done, but a few less. If there is need to update a number of points
it is thus more efficient to do so in a group update:
```Julia
setab(mytriangle, Point(1.7, 1.7), Point(1.3, 1.1))
```
combinations for all points exist. The name always contains the point names
in alphabetical order.

