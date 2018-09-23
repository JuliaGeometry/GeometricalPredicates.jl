# GeometricalPredicates.jl

[![Build Status](https://travis-ci.org/JuliaGeometry/GeometricalPredicates.jl.svg?branch=master)](https://travis-ci.org/JuliaGeometry/GeometricalPredicates.jl)
[![GeometricalPredicates](http://pkg.julialang.org/badges/GeometricalPredicates_0.6.svg)](http://pkg.julialang.org/detail/GeometricalPredicates)
[![Coverage Status](https://coveralls.io/repos/JuliaGeometry/GeometricalPredicates.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaGeometry/GeometricalPredicates.jl?branch=master)

Fast, robust 2D and 3D geometrical predicates on generic point types.
Implementation follows algorithms described in the [Arepo paper](http://arxiv.org/abs/0901.4107)
and used (for e.g.) in the [Illustris Simulation](http://www.illustris-project.org/). License: MIT. Bug reports welcome!

How does it work?
--------------------
Calculations are initially performed on `Float64` while bounding max
absolute errors. If unable to determine result, fall back to exact
calculation using `BigInt`s. This is a form of floating point filtering.
Most calculations are cached for fast repeated testing of
incircle/intriangle predicates.

Current limitations
--------------------
* Due to the numerical methods used, all coordinates are internally represented as `Float64`. In addition all must reside in the range `1.0<=x<2.0`. In this range, according to IEEE754, `Float64`s have a constant exponent, hence their mantissa can be used for a one to one mapping to integers, which in turn are used for exact calculations using `BigInt`s.
* It is assumed that primitive vertices don't overlap. It is currently the responsibility of the user to make sure this is the case.
* It is assumed tetrahedron vertices don't all lie in the same line. It is the user's responsibility to make sure it is so.
* Testing points are assumed not to overlap any vertices. As usual, it is the user's responsibility to make sure this is the case.
Except for the 1st restriction, all others could be easily implemented. Currently these features are not needed by me. If you need missing features, please open an issue I may develop it!

How to use?
--------------
### Installation
```julia
]add GeometricalPredicates
```
For Julia 0.6 and older
```julia
Pkg.add("GeometricalPredicates")
```

### Points
```julia
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
```julia
type MyCustomPointType <: AbstractPoint2D
    _x::FLoat64
    _y::Float64
    _mass::Float64
end

getx(p::MyCustomPointType) = p._x
gety(p::MyCustomPointType) = p._y
```
implementing `getx`, `gety`, and `getz` for 3D points is necessary
as this is the interface the package is expecting. These function should return `Float64`.
Points can be either immutables or types. Default `Point2D` and `Point3D` are immutables.

The point coordinates must reside in a region `1.0 <= x < 2.0`. Read above on
why this limitation is necessary. For convenience there are 2 constants defined,
`min_coord` and `max coord` representing the minimal and maximal feasible values
of coordinates.

### Triangles and Tetrahedrons (..aka Primitives)
A triangle is the 2D primitive, and a tetrahedron is the 3D primitive.
```julia
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
points within are updated. Read below for the definition of orientation. The
triangle could be created using any points inheriting from `AbstractPoint2D`
which implement `getx` and `gety`, or using coordinates directly:
```julia
mytriangle = Primitive(1.1, 1.1, 1.9, 1.1, 1.1, 1.9)

# Getting point `a` in the triangle
geta(mytriangle) # -> Point2D(1.1, 1.1)
getb(mytriangle) # -> Point2D(1.9, 1.1)
getc(mytriangle) # -> Point2D(1.1, 1.9)
```
The same goes for tetrahedrons, except we now use 4 3D points instead of 3 2D ones:
```julia
# create a tetrahedron using 4 points
a = Point(1.1, 1.1, 1.1)
b = Point(1.9, 1.1, 1.1)
c = Point(1.1, 1.9, 1.1)
d = Point(1.1, 1.1, 1.9)
mytetraedron = Primitive(a, b, c, d)
typeof(mytetrahedron) # -> UnOrientedTetrahedron{Point3D}
```
For certain applications we use primitives with known orientation. In those cases
there should be no need to calculate it. This is achieved by passing an `orientation`
flag to the `Primitive` creation function:
```julia
mytetrahedron = Primitive(a, b, c, d, positivelyoriented)
typeof(mytetrahedron) # -> PositivelyOrientedTetrahedron{Point3D}
orientation(mytetrahedron) # -> constant 1, not calculated
mytetrahedron = Primitive(a, b, c, d, negativelyoriented)
typeof(mytetrahedron) # -> NegativelyOrientedTetrahedron{Point3D}
orientation(mytetrahedron) # -> constant -1, not calculated
```
Note that when the primitive is oriented the real orientation is never calculated.
It is assumed that the user knows what he's doing. If in doubt, just use unoriented
primitives at the cost of actual calculation.

Updating points in primitives can be done with `seta`, `setb`, etc. methods:
```julia
seta(mytriangle, Point(1.7, 1.7))
```
Updating a point in a primitive will fire all relevant pre-calculations. i.e. if the triangle
is unoriented then orientation will be calculated. If it is oriented then still a few other
pre-calculations will be done, but a few less. If there is need to update a number of points
it is thus more efficient to do so in a group update:
```julia
setab(mytriangle, Point(1.7, 1.7), Point(1.3, 1.1))
setbcd(mytetrahedron, Point(1.1, 1.1, 1.2), Point(1.2,1.1,1.3), Point(1.4,1.1,1.2))
```
combinations for all points exist. The name always contains the point names
in alphabetical order. As long as inner primitive data is not changed manually, it will
keep giving correct results for all functions in this package.

### Predicates
`incircle` is the popular name to test whether a point lies inside of the sphere
defined by the primitive points:
```julia
a = Point(1.1, 1.1)
b = Point(1.5, 1.1)
c = Point(1.1, 1.5)
mytriangle = Primitive(a, b, c)
incircle(mytriangle, Point(1.9, 1.9)) # -> -1, i.e. outside
incircle(mytriangle, Point(1.2, 1.2)) # -> +1, i.e. inside
incircle(mytriangle, Point(1.5, 1.5)) # ->  0, i.e. point is exactly on circle
```
There is one more possibility. If the circle defined by our primitive has infinite radius
then it is impossible to tell whether the point is inside or outside:
```julia
a = Point(1.1, 1.1)
b = Point(1.2, 1.2)
c = Point(1.3, 1.3)
mytriangle = Primitive(a, b, c)
incircle(mytriangle, Point(1.3, 1.4)) # -> +2, i.e. cannot tell
```

`intriangle` is a popular name to test whether a point lies inside of the primitive:
```julia
a = Point(1.1, 1.1)
b = Point(1.5, 1.1)
c = Point(1.1, 1.5)
mytriangle = Primitive(a, b, c)
intriangle(mytriangle, Point(1.2, 1.2)) # -> +1, i.e. inside
intriangle(mytriangle, Point(1.6, 1.6)) # -> -1, i.e. outside
intriangle(mytriangle, Point(1.3, 1.1)) # ->  4, i.e. exactly on ab
intriangle(mytriangle, Point(1.1, 1.3)) # ->  3, i.e. exactly on ac
intriangle(mytriangle, Point(1.3, 1.3)) # ->  2, i.e. exactly on bc

```
Here any negative number means outside. The exact value gives some information regarding
the direction in which the point lies outside:
* `-1` means the test point is in front of a, and outside of the triangle
* `-2` means the test point is in front of b, and outside of the triangle
* `-4` means the test point is in front of c, and outside of the triangle
same goes for tetrahedrons. Note that the point could be both in front of a and b. In
cases as this one is arbitrarily chosen, all in name of performance.

`1` means test point is inside. But there are other possible positive values:
* `1 + 1 = 2` means the test point is in front of a, exactly on the triangle
* `1 + 2 = 3` means the test point is in front of b, exactly on the triangle
* `1 + 3 = 4` means the test point is in front of c, exactly on the triangle

same extends for tetrahedrons.

### Lines and Polygons

`orientation` tests for the primitive orientation. In 2D this means:
* ` 1` --> point `c` is to the left of the line defined by `ab` (with directionality from `a` to `b`)
* `-1` --> point `c` is to the right
* ` 0` --> point `c` is exactly on the line

in 3D it means:
* ` 1` --> point `d` is above the plane defined by `abc` (note "above" here means the direction of the plane normal, which depends on its orientation)
* `-1` --> point `d` is below the plane
* ` 0` --> point `c` is exactly on the plane

Another convenience API to test for orientation in 2D is using a line. It has some performance advantages over creating a triangle:
```julia
a = Point(1.1, 1.1)
b = Point(1.5, 1.5)

l = Line(a, b)
println(orientation(l, Point(1.4, 1.6))) # -->  1
println(orientation(l, Point(1.4,1.05))) # --> -1
println(orientation(l, Point(1.4,1.40))) # -->  0
```

One can also create simple 2D polygons:
```julia
ll = Point(1.0,1.0)
lr = Point(1.2,1.0)
ur = Point(1.2,1.2)
ul = Point(1.0,1.2)
poly = Polygon(ll, lr, ur, ul)
inpolygon(poly, Point(1.1,1.1))  # assumes it is convex
```


### Basic geometrical properties
`orientation` gives the primitive orientation. `area`, `volume`, `centroid`, `circumcenter`, `circumradius2` are all exported and I hope self descriptive.

### Spatial ordering
Scale and scale-free Peano-Hilbert ordering is available. Look [here](http://doc.cgal.org/latest/Spatial_sorting/index.html) for a nice explanation on Hilbert sorting and [here](http://doc.cgal.org/latest/Spatial_sorting/classCGAL_1_1Multiscale__sort.html) for a nice explanation of multiscale sort. Both are implemented here:

```julia
p = Point(1.1, 1.2)
peanokey(p, 4) # -> 6, the peano key in a regular grid of 2^4 X 2^4 cells

p = Point(1.1, 1.2, 1.3)
peanokey(p, 4) # -> 94, the peano key in a regular grid of 2^4 X 2^4 X 2^4 cells
```

The number of cells doesn't need to be specified. The default for 2D is `2^31 X 2^31` and for 3D `2^21 X 2^21 X 2^21`.
You can also do the inverse, and get a point from a key:
```julia
Point2D(6, 4) # -> Point2D(1.0625,1.1875)
```
in a finer grid we would get back something more accurate.

So that was scale-dependent Hilbert stuff, which is good to say balance workload between computing nodes.
When you need to order particles spatially it is better to use a scale independent method, like the Hilbert ordering here:

```julia
a = [Point(1.0+rand(), 1.0+rand() for i in 1:1e6]
hilbertsort!(a)
```

Here keys are never calculated, and there is no grid, it uses a median policy, adjusting to the actual
point distribution. There are a few parameters with good defaults, see links above to understand what they mean.
For an algorithm such a Delaunay tesselation it is sometimes better to use a multi-scale sort with a Hilbert sort, like this:
```julia
mssort!(a)
```
of course this adds a few more default parameters, again with decent defaults, read above links to understand.
