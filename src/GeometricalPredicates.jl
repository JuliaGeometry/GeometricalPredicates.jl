VERSION >= v"0.4.0-dev+6521" && __precompile__()

module GeometricalPredicates

# Fast, robust 2D and 3D geometrical predicates on generic point types.
# Implementation follows algorithms described in http://arxiv.org/abs/0901.4107
# and used (for e.g.) in the Illustris Simulation
# http://www.illustris-project.org/

# Author: Ariel Keselman (skariel@gmail.com)
# License: MIT
# Bug reports welcome!

using Compat

export
    min_coord, max_coord,

    AbstractPoint,
        AbstractPoint2D,
        AbstractPoint3D,

    AbstractPositivelyOrientedPrimitive,
        AbstractPositivelyOrientedTriangle,
        AbstractPositivelyOrientedTetrahedron,
    AbstractNegativelyOrientedPrimitive,
        AbstractNegativelyOrientedTriangle,
        AbstractNegativelyOrientedTetrahedron,
    AbstractTriangleUnOriented,
    AbstractTetrahedronUnOriented,

    TriangleTypes, TetrahedronTypes,

    positivelyoriented, negativelyoriented, unoriented, orientation,

    Point, Point2D, Point3D, Line, Polygon, getx, gety, getz, geta, getb, getc, getd,
    seta, setb, setc, setd, setabc, setabcd,
    setab, setbc, setcd, setac, setad, setbd,
    setabd, setacd, setbcd,
    getlines, getpoints,

    Line2D, Polygon2D, Primitive, Triangle, Tetrahedron,

    length2, area, volume, centroid, circumcenter, circumradius2, incircle, intriangle, inpolygon,

    peanokey, hilbertsort!, mssort!, clean!


const _float_err = eps(Float64)
const _abs_err_incircle_2d = 12*_float_err
const _abs_err_incircle_3d = 48*_float_err
const _abs_err_orientation_2d = 2*_float_err
const _abs_err_orientation_3d = 6*_float_err
const _abs_err_intriangle = 6*_float_err
const _abs_err_intriangle_zero = 2*_float_err
const _abs_err_intetra = 24*_float_err
const _abs_err_intetra_zero = 6*_float_err

const min_coord = 1.0
const max_coord = 2.0 - eps(Float64)

abstract AbstractPoint
abstract AbstractPoint2D <: AbstractPoint
abstract AbstractPoint3D <: AbstractPoint

abstract AbstractLine2D
abstract AbstractPolygon2D

abstract AbstractPrimitive
abstract AbstractUnOrientedPrimitive <: AbstractPrimitive
abstract AbstractOrientedPrimitive <: AbstractPrimitive
abstract AbstractPositivelyOrientedPrimitive <: AbstractOrientedPrimitive
abstract AbstractNegativelyOrientedPrimitive <: AbstractOrientedPrimitive

abstract AbstractTriangleUnOriented <: AbstractUnOrientedPrimitive
abstract AbstractTetrahedronUnOriented <: AbstractUnOrientedPrimitive
abstract AbstractPositivelyOrientedTriangle <: AbstractPositivelyOrientedPrimitive
abstract AbstractNegativelyOrientedTriangle <: AbstractNegativelyOrientedPrimitive
abstract AbstractPositivelyOrientedTetrahedron <: AbstractPositivelyOrientedPrimitive
abstract AbstractNegativelyOrientedTetrahedron <: AbstractNegativelyOrientedPrimitive

typealias TriangleTypes Union{AbstractTriangleUnOriented, AbstractPositivelyOrientedTriangle, AbstractNegativelyOrientedTriangle}
typealias TetrahedronTypes Union{AbstractTetrahedronUnOriented, AbstractPositivelyOrientedTetrahedron, AbstractNegativelyOrientedTetrahedron}

# standard 2D point
immutable Point2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    Point2D(x::Float64,y::Float64) = new(x, y)
end
Point2D() = Point2D(0., 0.)

getx(p::Point2D) = p._x
gety(p::Point2D) = p._y

# standard 3D point
immutable Point3D <: AbstractPoint3D
    _x::Float64
    _y::Float64
    _z::Float64
    Point3D(x::Float64,y::Float64,z::Float64) = new(x, y, z)
end
Point3D() = Point3D(0., 0., 0.)

getx(p::Point3D) = p._x
gety(p::Point3D) = p._y
getz(p::Point3D) = p._z

Point(x::Real, y::Real) = Point2D(x, y)
Point(x::Real, y::Real, z::Real) = Point3D(x, y, z)

immutable Line2D{T<:AbstractPoint2D} <: AbstractLine2D
    _a::T
    _b::T
    _bx::Float64
    _by::Float64
    function Line2D(a::T, b::T)
        const bx = getx(b) - getx(a)
        const by = gety(b) - gety(a)
        new(a, b, bx, by)
    end
end

Line2D{T<:AbstractPoint2D}(a::T, b::T) = Line2D{T}(a, b)

Line{T<:AbstractPoint2D}(a::T, b::T) = Line2D(a, b)

geta(l::Line2D) = l._a
getb(l::Line2D) = l._b

length2(l::Line2D) = l._bx*l._bx + l._by*l._by


# fine filtered orientation
function _sz_orientation(l::Line2D, p::AbstractPoint2D)
    const cx = getx(p) - getx(geta(l))
    const cy = gety(p) - gety(geta(l))
    const _pr2 = -l._bx*cy + l._by*cx
    const sz = abs(cx) + abs(cy) + abs(l._bx) + abs(l._by)
    if _pr2 < -_abs_err_orientation_2d*sz
        1
    elseif _pr2 > _abs_err_orientation_2d*sz
        -1
    else
        _exact_sign_orientation_determinant!(
            _extract_bigint(getx(geta(l))), _extract_bigint(gety(geta(l))),
            _extract_bigint(getx(getb(l))), _extract_bigint(gety(getb(l))),
            _extract_bigint(getx(p)), _extract_bigint(gety(p)))
    end
end

# gross filtered orientation, asumming maximal line size (=1.0)
function orientation(l::Line2D, p::AbstractPoint2D)
    const cx = getx(p) - getx(geta(l))
    const cy = gety(p) - gety(geta(l))
    const _pr2 = -l._bx*cy + l._by*cx
    if _pr2 < -_abs_err_orientation_2d
        1
    elseif _pr2 > _abs_err_orientation_2d
        -1
    else
        _sz_orientation(l, p)
    end
end

"a simple polygon"
immutable Polygon2D{T<:AbstractPoint2D} <: AbstractPolygon2D
    _p::Vector{T}
    _l::Vector{AbstractLine2D}
    function Polygon2D(p::T...)
        l=AbstractLine2D[]
        for i=1:length(p)-1
            push!(l,Line(p[i],p[i+1]))
        end
        push!(l,Line(p[end],p[1]))
        new([p...],l)
    end
end

Polygon2D{T<:AbstractPoint2D}(p::T...) = Polygon2D{T}(p...)

Polygon{T<:AbstractPoint2D}(p::T...) = Polygon2D(p...)

"return the points of a Polygon"
getpoints(polygon::Polygon2D) = polygon._p

"return the lines of a Polygon"
getlines(polygon::Polygon2D) = polygon._l

"return true if the Point is inside the Polygon, which is assumed to be convex"
function inpolygon(polygon::Polygon2D, point::AbstractPoint2D)
    lines = getlines(polygon)
    side = orientation(lines[1], point)
    for i = 2:length(lines)
        orientation(lines[i], point) == side || return false
    end
    true
end

macro _define_triangle_type(name, abstracttype)
    oriented = !contains(string(name), "UnOriented")
    esc(parse("""
        type $name{T<:AbstractPoint2D} <: $abstracttype
            _a::T; _b::T; _c::T
            _bx::Float64; _by::Float64
            _cx::Float64; _cy::Float64
            _px::Float64; _py::Float64
            _pr2::Float64
            $(oriented ? "" : "_o::Int8")
            function $name(a::T, b::T, c::T)
                t = new(a, b, c, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0$(oriented? "":", 0"))
                clean!(t)
                t
            end
        end
    """))
end

@_define_triangle_type(UnOrientedTriangle, AbstractTriangleUnOriented)
@_define_triangle_type(PositivelyOrientedTriangle, AbstractPositivelyOrientedTriangle)
@_define_triangle_type(NegativelyOrientedTriangle, AbstractNegativelyOrientedTriangle)

abstract AbstractOrientation
type PositivelyOriented <: AbstractOrientation; end
type NegativelyOriented <: AbstractOrientation; end
type UnOriented <: AbstractOrientation; end

const positivelyoriented = PositivelyOriented()
const negativelyoriented = NegativelyOriented()
const unoriented = UnOriented()

Triangle{T<:AbstractPoint2D}(a::T, b::T, c::T, ::PositivelyOriented) = PositivelyOrientedTriangle{T}(a, b, c)
Triangle{T<:AbstractPoint2D}(a::T, b::T, c::T, ::NegativelyOriented) = NegativelyOrientedTriangle{T}(a, b, c)
Triangle{T<:AbstractPoint2D}(a::T, b::T, c::T, ::UnOriented) = UnOrientedTriangle{T}(a, b, c)
Triangle{T<:AbstractPoint2D}(a::T, b::T, c::T) = Triangle(a, b, c, unoriented)
Triangle(ax::Float64, ay::Float64, bx::Float64, by::Float64, cx::Float64, cy::Float64, orientation::AbstractOrientation=unoriented) =
    Triangle(Point2D(ax, ay), Point2D(bx, by), Point2D(cx, cy), orientation)
Primitive(ax::Float64, ay::Float64, bx::Float64, by::Float64, cx::Float64, cy::Float64, orientation::AbstractOrientation=unoriented) =
    Triangle(Point2D(ax, ay), Point2D(bx, by), Point2D(cx, cy), orientation)
Primitive{T<:AbstractPoint2D}(a::T, b::T, c::T, orientation::AbstractOrientation=unoriented) =
    Triangle(a, b, c, orientation)


area(tr::TriangleTypes) = abs(tr._pr2)/2

centroid(tr::TriangleTypes) =
    Point2D(
        (getx(geta(tr)) + getx(getb(tr)) + getx(getc(tr))) / 3.0,
        (gety(geta(tr)) + gety(getb(tr)) + gety(getc(tr))) / 3.0)

function circumcenter(tr::TriangleTypes)
    const d = -2.0 * tr._pr2
    Point2D(
        tr._px/d + getx(geta(tr)),
        tr._py/d + gety(geta(tr))
        )
end

function circumradius2(tr::TriangleTypes)
    const c = circumcenter(tr)
    const x = getx(c) - getx(geta(tr))
    const y = gety(c) - gety(geta(tr))
    x*x + y*y
end

macro _define_tetrahedron_type(name, abstracttype)
    oriented = !contains(string(name), "UnOriented")
    esc(parse("""
        type $name{T<:AbstractPoint3D} <: $abstracttype
            _a::T; _b::T; _c::T; _d::T
            _bx::Float64; _by::Float64; _bz::Float64
            _cx::Float64; _cy::Float64; _cz::Float64
            _dx::Float64; _dy::Float64; _dz::Float64
            _px::Float64; _py::Float64; _pz::Float64
            _pr2::Float64
            $(oriented? "":"_o::Int8")
            function $name(a::AbstractPoint3D, b::AbstractPoint3D, c::AbstractPoint3D, d::AbstractPoint3D)
                t = new(a, b, c, d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0$(oriented? "" : ", 0"))
                clean!(t)
                t
            end
        end
    """))
end

@_define_tetrahedron_type(UnOrientedTetrahedron, AbstractTetrahedronUnOriented)
@_define_tetrahedron_type(PositivelyOrientedTetrahedron, AbstractPositivelyOrientedTetrahedron)
@_define_tetrahedron_type(NegativelyOrientedTetrahedron, AbstractNegativelyOrientedTetrahedron)

Tetrahedron{T<:AbstractPoint3D}(a::T, b::T, c::T, d::T, ::PositivelyOriented) = PositivelyOrientedTetrahedron{T}(a, b, c, d)
Tetrahedron{T<:AbstractPoint3D}(a::T, b::T, c::T, d::T, ::NegativelyOriented) = NegativelyOrientedTetrahedron{T}(a, b, c, d)
Tetrahedron{T<:AbstractPoint3D}(a::T, b::T, c::T, d::T, ::UnOriented) = UnOrientedTetrahedron{T}(a, b, c, d)
Tetrahedron{T<:AbstractPoint3D}(a::T, b::T, c::T, d::T) = Tetrahedron(a, b, c, d, unoriented)
Tetrahedron(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64,
    cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64, orientation::AbstractOrientation=unoriented) =
        Tetrahedron(Point3D(ax,ay,az), Point3D(bx,by,bz), Point3D(cx,cy,cz), Point3D(dx,dy,dz), orientation)
Primitive(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64,
    cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64, orientation::AbstractOrientation=unoriented) =
        Tetrahedron(Point3D(ax,ay,az), Point3D(bx,by,bz), Point3D(cx,cy,cz), Point3D(dx,dy,dz), orientation)
Primitive{T<:AbstractPoint3D}(a::T, b::T, c::T, d::T, orientation::AbstractOrientation=unoriented) =
        Tetrahedron(a, b, c, d, orientation)

volume(tr::TetrahedronTypes) = abs(tr._pr2)/2

centroid(tr::TetrahedronTypes) =
    Point3D(
        (getx(geta(tr)) + getx(getb(tr)) + getx(getc(tr)) + getx(getd(tr))) / 4.0,
        (gety(geta(tr)) + gety(getb(tr)) + gety(getc(tr)) + gety(getd(tr))) / 4.0,
        (getz(geta(tr)) + getz(getb(tr)) + getz(getc(tr)) + getz(getd(tr))) / 4.0)

function circumcenter(tr::TetrahedronTypes)
    const d = -2.0 * tr._pr2
    Point3D(
        tr._px/d + getx(geta(tr)),
        tr._py/d + gety(geta(tr)),
        tr._pz/d + getz(geta(tr))
        )
end

function circumradius2(tr::TetrahedronTypes)
    const c = circumcenter(tr)
    const x = getx(c) - getx(geta(tr))
    const y = gety(c) - gety(geta(tr))
    const z = getz(c) - getz(geta(tr))
    x*x + y*y + z*z
end

# extract exact integer representation of float to be used in exact calculations when needed
_extract_int(n::Float64) = reinterpret(UInt64, n) & 0x000fffffffffffff
_extract_bigint(n::Float64) = BigInt(_extract_int(n))

# functions to re-validate cached pre-calculations in primitives
function _clean!(t::TriangleTypes)
    t._bx = getx(getb(t))-getx(geta(t)); t._by = gety(getb(t))-gety(geta(t));
    t._cx = getx(getc(t))-getx(geta(t)); t._cy = gety(getc(t))-gety(geta(t));

    const br2 = t._bx*t._bx+t._by*t._by
    const cr2 = t._cx*t._cx+t._cy*t._cy

    t._px  =  br2*t._cy - t._by*cr2
    t._py  = -br2*t._cx + t._bx*cr2
    t._pr2 = -t._bx *t._cy + t._by*t._cx
end

function _sz_clean!(t::AbstractTriangleUnOriented)
    sz = abs(t._bx)+abs(t._by)+abs(t._cx)+abs(t._cy)
    if t._pr2 < -_abs_err_orientation_2d*sz
        t._o = 1
    elseif t._pr2 > _abs_err_orientation_2d*sz
        t._o = -1
    else
        t._o = _exact_sign_orientation_determinant!(
            _extract_bigint(getx(geta(t))), _extract_bigint(gety(geta(t))),
            _extract_bigint(getx(getb(t))), _extract_bigint(gety(getb(t))),
            _extract_bigint(getx(getc(t))), _extract_bigint(gety(getc(t))))
    end
end
function clean!(t::AbstractTriangleUnOriented)
    _clean!(t)
    if t._pr2 < -_abs_err_orientation_2d
        t._o = 1
    elseif t._pr2 > _abs_err_orientation_2d
        t._o = -1
    else
        _sz_clean!(t)
    end
end

function _clean!(t::TetrahedronTypes)
    t._bx = getx(getb(t))-getx(geta(t)); t._by = gety(getb(t))-gety(geta(t)); t._bz = getz(getb(t))-getz(geta(t))
    t._cx = getx(getc(t))-getx(geta(t)); t._cy = gety(getc(t))-gety(geta(t)); t._cz = getz(getc(t))-getz(geta(t))
    t._dx = getx(getd(t))-getx(geta(t)); t._dy = gety(getd(t))-gety(geta(t)); t._dz = getz(getd(t))-getz(geta(t))

    const br2 = t._bx*t._bx+t._by*t._by+t._bz*t._bz
    const cr2 = t._cx*t._cx+t._cy*t._cy+t._cz*t._cz
    const dr2 = t._dx*t._dx+t._dy*t._dy+t._dz*t._dz

    t._px  =  t._by*t._cz*dr2 - br2*t._cz*t._dy - t._by*cr2*t._dz - t._bz*t._cy*dr2 + br2*t._cy*t._dz + t._bz*cr2*t._dy
    t._py  =  br2*t._cz*t._dx + t._bz*t._cx*dr2 - br2*t._cx*t._dz - t._bx*t._cz*dr2 - t._bz*cr2*t._dx + t._bx*cr2*t._dz
    t._pz  =  br2*t._cx*t._dy + t._bx*t._cy*dr2 + t._by*cr2*t._dx - br2*t._cy*t._dx - t._bx*cr2*t._dy - t._by*t._cx*dr2
    t._pr2 = -t._bx*t._cy*t._dz  + t._bx*t._cz*t._dy  + t._by*t._cx*t._dz  - t._by*t._cz*t._dx  - t._bz*t._cx*t._dy  + t._bz*t._cy*t._dx
end

function _sz_clean!(t::AbstractTetrahedronUnOriented)
    const sz = abs(t._bx) + abs(t._by) + abs(t._bz) +
                    abs(t._cx) + abs(t._cy) + abs(t._cz) +
                    abs(t._dx) + abs(t._dy) + abs(t._dz)

    # calculate the orientation
    if t._pr2 < -_abs_err_orientation_3d*sz
        t._o = 1
    elseif t._pr2 > _abs_err_orientation_3d*sz
        t._o = -1
    else
        # exact calculation is required
        t._o = _exact_sign_orientation_determinant!(
            _extract_bigint(getx(geta(t))), _extract_bigint(gety(geta(t))), _extract_bigint(getz(geta(t))),
            _extract_bigint(getx(getb(t))), _extract_bigint(gety(getb(t))), _extract_bigint(getz(getb(t))),
            _extract_bigint(getx(getc(t))), _extract_bigint(gety(getc(t))), _extract_bigint(getz(getc(t))),
            _extract_bigint(getx(getd(t))), _extract_bigint(gety(getd(t))), _extract_bigint(getz(getd(t))))
    end
end

function clean!(t::AbstractTetrahedronUnOriented)
    _clean!(t)
    # calculate the orientation
    if t._pr2 < -_abs_err_orientation_3d
        t._o = 1
    elseif t._pr2 > _abs_err_orientation_3d
        t._o = -1
    else
        _sz_clean!(t)
    end
end

clean!(t::AbstractOrientedPrimitive) = _clean!(t)

# getting points of the primitive
geta(t::AbstractPrimitive) = t._a
getb(t::AbstractPrimitive) = t._b
getc(t::AbstractPrimitive) = t._c
getd(t::TetrahedronTypes) = t._d

# changing points in the primitive
seta(t::AbstractPrimitive, p::AbstractPoint) = (t._a=p; clean!(t))
setb(t::AbstractPrimitive, p::AbstractPoint) = (t._b=p; clean!(t))
setc(t::AbstractPrimitive, p::AbstractPoint) = (t._c=p; clean!(t))
setd(t::TetrahedronTypes, p::AbstractPoint3D) = (t._d=p; clean!(t))
setabc(t::AbstractPrimitive, pa::AbstractPoint, pb::AbstractPoint, pc::AbstractPoint) = (t._a=pa; t._b=pb; t._c=pc; clean!(t))

setab(t::AbstractPrimitive, pa::AbstractPoint, pb::AbstractPoint) = (t._a=pa; t._b=pb; clean!(t))
setbc(t::AbstractPrimitive, pb::AbstractPoint, pc::AbstractPoint) = (t._b=pb; t._c=pc; clean!(t))
setac(t::AbstractPrimitive, pa::AbstractPoint, pc::AbstractPoint) = (t._a=pa; t._c=pc; clean!(t))

setabcd(t::TetrahedronTypes, pa::AbstractPoint3D, pb::AbstractPoint3D, pc::AbstractPoint3D, pd::AbstractPoint3D) = (t._a=pa; t._b=pb; t._c=pc; t._d=pd; clean!(t))

setabc(t::TetrahedronTypes, pa::AbstractPoint3D, pb::AbstractPoint3D, pc::AbstractPoint3D) = (t._a=pa; t._b=pb; t._c=pc; clean!(t))
setabd(t::TetrahedronTypes, pa::AbstractPoint3D, pb::AbstractPoint3D, pd::AbstractPoint3D) = (t._a=pa; t._b=pb; t._d=pd; clean!(t))
setacd(t::TetrahedronTypes, pa::AbstractPoint3D, pc::AbstractPoint3D, pd::AbstractPoint3D) = (t._a=pa; t._c=pc; t._d=pd; clean!(t))
setbcd(t::TetrahedronTypes, pb::AbstractPoint3D, pc::AbstractPoint3D, pd::AbstractPoint3D) = (t._b=pb; t._c=pc; t._d=pd; clean!(t))

setad(t::TetrahedronTypes, pa::AbstractPoint3D, pd::AbstractPoint3D) = (t._a=pa; t._d=pd; clean!(t))
setbd(t::TetrahedronTypes, pb::AbstractPoint3D, pd::AbstractPoint3D) = (t._b=pb; t._d=pd; clean!(t))
setcd(t::TetrahedronTypes, pc::AbstractPoint3D, pd::AbstractPoint3D) = (t._c=pc; t._d=pd; clean!(t))

orientation(p::AbstractUnOrientedPrimitive) = p._o
orientation(p::AbstractPositivelyOrientedPrimitive) = 1
orientation(p::AbstractNegativelyOrientedPrimitive) = -1
orientation(ax::Float64, ay::Float64, bx::Float64, by::Float64, cx::Float64, cy::Float64) =
    orientation(Triangle(ax, ay, bx, by, cx, cy))
orientation(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64, cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64) =
    orientation(Tetrahedron(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz))

# exact orientation for triangle
function _exact_sign_orientation_determinant!(ax::BigInt, ay::BigInt, bx::BigInt, by::BigInt, cx::BigInt, cy::BigInt)
    bx -= ax; by -= ay
    cx -= ax; cy -= ay
    @compat Int64(sign(bx*cy - by*cx))
end

# exact orientation for tetrahedron
function _exact_sign_orientation_determinant!(ax::BigInt, ay::BigInt, az::BigInt, bx::BigInt, by::BigInt, bz::BigInt, cx::BigInt, cy::BigInt, cz::BigInt, dx::BigInt, dy::BigInt, dz::BigInt)
    bx -= ax; by -= ay; bz -= az
    cx -= ax; cy -= ay; cz -= az
    dx -= ax; dy -= ay; dz -= az
    @compat Int64(sign(+bx*cy*dz - bx*cz*dy - by*cx*dz + by*cz*dx + bz*cx*dy - bz*cy*dx))
end

# exact incircle for triangle
function _exact_sign_incircle_determinant!(ax::BigInt, ay::BigInt, bx::BigInt, by::BigInt, cx::BigInt, cy::BigInt, px::BigInt, py::BigInt)
    bx -= ax; by -= ay;
    cx -= ax; cy -= ay;
    px -= ax; py -= ay;
    const br2 = bx*bx+by*by
    const cr2 = cx*cx+cy*cy
    const pr2 = px*px+py*py
    @compat Int64(sign(-br2*cx*py + br2*cy*px + bx*cr2*py - bx*cy*pr2 - by*cr2*px + by*cx*pr2))
end

# exact incircle for tetrahedron
function _exact_sign_incircle_determinant!(ax::BigInt, ay::BigInt, az::BigInt, bx::BigInt, by::BigInt, bz::BigInt, cx::BigInt, cy::BigInt, cz::BigInt, dx::BigInt, dy::BigInt, dz::BigInt, px::BigInt, py::BigInt, pz::BigInt)
    bx -= ax; by -= ay; bz -= az
    cx -= ax; cy -= ay; cz -= az
    dx -= ax; dy -= ay; dz -= az
    px -= ax; py -= ay; pz -= az
    const br2 = bx*bx+by*by+bz*bz
    const cr2 = cx*cx+cy*cy+cz*cz
    const dr2 = dx*dx+dy*dy+dz*dz
    const pr2 = px*px+py*py+pz*pz
    @compat Int64(sign(
        +br2*cx*dy*pz - br2*cx*dz*py - br2*cy*dx*pz + br2*cy*dz*px +
         br2*cz*dx*py - br2*cz*dy*px - bx*cr2*dy*pz + bx*cr2*dz*py +
         bx*cy*dr2*pz - bx*cy*dz*pr2 - bx*cz*dr2*py + bx*cz*dy*pr2 +
         by*cr2*dx*pz - by*cr2*dz*px - by*cx*dr2*pz + by*cx*dz*pr2 +
         by*cz*dr2*px - by*cz*dx*pr2 - bz*cr2*dx*py + bz*cr2*dy*px +
         bz*cx*dr2*py - bz*cx*dy*pr2 - bz*cy*dr2*px + bz*cy*dx*pr2))
end

# finer filtered incircle
function _sz_incircle(t::TriangleTypes, p::AbstractPoint2D, px::Float64, py::Float64, pr2::Float64)
    if orientation(t) != 0
        sz = abs(px)+abs(py)+abs(t._bx)+abs(t._by)+abs(t._cx)+abs(t._cy)
        d = t._px*px + t._py*py + t._pr2*pr2
        if d < -_abs_err_incircle_2d*sz
            return -orientation(t)
        elseif d > _abs_err_incircle_2d*sz
            return orientation(t)
        end
    end

    const exact_in = _exact_sign_incircle_determinant!(
        _extract_bigint(getx(geta(t))), _extract_bigint(gety(geta(t))),
        _extract_bigint(getx(getb(t))), _extract_bigint(gety(getb(t))),
        _extract_bigint(getx(getc(t))), _extract_bigint(gety(getc(t))),
        _extract_bigint(getx(p))  , _extract_bigint(gety(p)))


    if orientation(t) != 0
        return orientation(t)*exact_in
    elseif exact_in == 0
        return 1
    else
        return 2
    end
end

# gross filtered incircle, asumming maximal triangle size (=1.0)
function incircle(t::TriangleTypes, p::AbstractPoint2D)
    px  = getx(p) - getx(geta(t))
    py  = gety(p) - gety(geta(t))
    pr2 = px*px + py*py
    if orientation(t) != 0
        d = t._px*px + t._py*py + t._pr2*pr2
        if d < -_abs_err_incircle_2d
            return -orientation(t)
        elseif d > _abs_err_incircle_2d
            return orientation(t)
        end
    end
    _sz_incircle(t, p, px, py, pr2)
end


# finer filtered incircle
function _sz_incircle(t::TetrahedronTypes, p::AbstractPoint3D)
    if orientation(t) != 0
        px  = getx(p) - getx(geta(t))
        py  = gety(p) - gety(geta(t))
        pz  = getz(p) - getz(geta(t))

        const sz = abs(px) + abs(py) + abs(pz) +
                        abs(t._bx) + abs(t._by) + abs(t._bz) +
                        abs(t._cx) + abs(t._cy) + abs(t._cz) +
                        abs(t._dx) + abs(t._dy) + abs(t._dz)

        pr2 = px*px + py*py + pz*pz
        d = t._px*px + t._py*py + t._pz*pz + t._pr2*pr2
        if d < -_abs_err_incircle_3d*sz
            return -orientation(t)
        elseif d > _abs_err_incircle_3d*sz
            return orientation(t)
        end
    end

    const exact_in = _exact_sign_incircle_determinant!(
        _extract_bigint(getx(geta(t))), _extract_bigint(gety(geta(t))), _extract_bigint(getz(geta(t))),
        _extract_bigint(getx(getb(t))), _extract_bigint(gety(getb(t))), _extract_bigint(getz(getb(t))),
        _extract_bigint(getx(getc(t))), _extract_bigint(gety(getc(t))), _extract_bigint(getz(getc(t))),
        _extract_bigint(getx(getd(t))), _extract_bigint(gety(getd(t))), _extract_bigint(getz(getd(t))),
        _extract_bigint(getx(p))  , _extract_bigint(gety(p))  , _extract_bigint(getz(p)))

    if orientation(t) != 0
        return orientation(t)*exact_in
    elseif exact_in == 0
        return 1
    else
        return 2
    end
end


# gross filtered incircle, asumming maximal possible tetra size (=1.0)
function incircle(t::TetrahedronTypes, p::AbstractPoint3D)
    if orientation(t) != 0
        px  = getx(p) - getx(geta(t))
        py  = gety(p) - gety(geta(t))
        pz  = getz(p) - getz(geta(t))
        pr2 = px*px + py*py + pz*pz
        d = t._px*px + t._py*py + t._pz*pz + t._pr2*pr2
        if d < -_abs_err_incircle_3d
            return -orientation(t)
        elseif d > _abs_err_incircle_3d
            return orientation(t)
        end
    end
    _sz_incircle(t, p)
end

# helper methods to use incircle directly with coordinates
incircle(ax::Float64, ay::Float64, bx::Float64, by::Float64, cx::Float64, cy::Float64, dx::Float64, dy::Float64) =
    incircle(Triangle(ax, ay, bx, by, cx, cy), Point2D(dx, dy))
incircle(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64, cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64, ex::Float64, ey::Float64, ez::Float64) =
    incircle(Tetrahedron(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz), Point3D(ex, ey, ez))


# exact intriangle (slow!)
function _exact_intriangle!(ax::BigInt, ay::BigInt, bx::BigInt, by::BigInt, cx::BigInt, cy::BigInt, px::BigInt, py::BigInt)
    bx -= ax; by -= ay;
    cx -= ax; cy -= ay;
    px -= ax; py -= ay;

    const nb = -cx*py + cy*px
    const nc = bx*py - by*px
    const denom = bx*cy - by*cx

    const sdenom = @compat Int64(sign(denom))
    if @compat Int64(sign(nb)) * sdenom < 0
        return -2
    end
    if @compat Int64(sign(nc)) * sdenom < 0
        return -3
    end
    const l = nb+nc - denom
    const sl = @compat Int64(sign(l)) * sdenom
    if sl > 0
        return -1
    end

    if @compat Int64(sign(nb)) == 0
        return 3
    end
    if @compat Int64(sign(nc)) == 0
        return 4
    end
    if sl == 0
        return 2
    end
    -sl
end

function _exact_intriangle!(ax::BigInt, ay::BigInt, az::BigInt, bx::BigInt, by::BigInt, bz::BigInt, cx::BigInt, cy::BigInt, cz::BigInt, dx::BigInt, dy::BigInt, dz::BigInt, px::BigInt, py::BigInt, pz::BigInt)
    bx -= ax; by -= ay; bz -= az
    cx -= ax; cy -= ay; cz -= az
    dx -= ax; dy -= ay; dz -= az
    px -= ax; py -= ay; pz -= az

    const denom = bx*cy*dz-bx*cz*dy-by*cx*dz+by*cz*dx+bz*cx*dy-bz*cy*dx

    const nb = cx*dy*pz-cx*dz*py-cy*dx*pz+cy*dz*px+cz*dx*py-cz*dy*px
    const sdenom = @compat Int64(sign(denom))
    if @compat Int64(sign(nb)) * sdenom < 0
        return -2
    end

    const nc = -bx*dy*pz+bx*dz*py+by*dx*pz-by*dz*px-bz*dx*py+bz*dy*px
    if @compat Int64(sign(nc)) * sdenom < 0
        return -3
    end

    const nd = bx*cy*pz-bx*cz*py-by*cx*pz+by*cz*px+bz*cx*py-bz*cy*px
    if @compat Int64(sign(nd)) * sdenom < 0
        return -4
    end

    const l = (nb+nc+nd - denom) * sdenom
    const sl = @compat Int64(sign(l))
    if sl > 0
        return -1
    end

    if @compat Int64(sign(nb)) == 0
        return 3
    end
    if @compat Int64(sign(nc)) == 0
        return 4
    end
    if @compat Int64(sign(nd)) == 0
        return 5
    end
    if sl == 0
        return 2
    end
    -sl
end


# helper methods to use the exact intriangle directly with coordinates
_exact_intriangle(t::TriangleTypes, p::AbstractPoint2D) =
    _exact_intriangle!(
        _extract_bigint(getx(geta(t))), _extract_bigint(gety(geta(t))),
        _extract_bigint(getx(getb(t))), _extract_bigint(gety(getb(t))),
        _extract_bigint(getx(getc(t))), _extract_bigint(gety(getc(t))),
        _extract_bigint(getx(p))  , _extract_bigint(gety(p)))

_exact_intriangle(t::TetrahedronTypes, p::AbstractPoint3D) =
    _exact_intriangle!(
        _extract_bigint(getx(geta(t))), _extract_bigint(gety(geta(t))), _extract_bigint(getz(geta(t))),
        _extract_bigint(getx(getb(t))), _extract_bigint(gety(getb(t))), _extract_bigint(getz(getb(t))),
        _extract_bigint(getx(getc(t))), _extract_bigint(gety(getc(t))), _extract_bigint(getz(getc(t))),
        _extract_bigint(getx(getd(t))), _extract_bigint(gety(getd(t))), _extract_bigint(getz(getd(t))),
        _extract_bigint(getx(p)), _extract_bigint(gety(p)), _extract_bigint(getz(p)))

# finer filter for intriangle, using triangle actual size
function _sz_intriangle(t::TriangleTypes, p::AbstractPoint2D, px::Float64, py::Float64)
    sz = abs(px)+abs(py)+abs(t._bx)+abs(t._by)+abs(t._cx)+abs(t._cy)

    const nb = (-t._cx*py + t._cy*px) * sign(-t._pr2)
    if nb < -_abs_err_intriangle_zero*sz
        return -2
    elseif nb < _abs_err_intriangle_zero*sz
        # we need an exact calculation
        return _exact_intriangle(t, p)
    end

    const nc = (t._bx*py - t._by*px) * sign(-t._pr2)
    if nc < -_abs_err_intriangle_zero*sz
        return -3
    elseif nc < _abs_err_intriangle_zero*sz
        # we need an exact calculation
        return _exact_intriangle(t, p)
    end

    const l = nb+nc + t._pr2*sign(-t._pr2)
    if l < -_abs_err_intriangle*sz
        return 1
    elseif l > _abs_err_intriangle*sz
        return -1
    end
    # we need an exact calculation
    _exact_intriangle(t, p)
end

# gross filter, assuming maximal possible triangle size (=1.0)
function intriangle(t::TriangleTypes, p::AbstractPoint2D)
    const px = getx(p) - getx(geta(t)); const py = gety(p) - gety(geta(t))

    const nb = (-t._cx*py + t._cy*px) * sign(-t._pr2)
    if nb < -_abs_err_intriangle_zero
        return -2
    elseif nb < _abs_err_intriangle_zero
        return _sz_intriangle(t, p, px, py)
    end

    const nc = (t._bx*py - t._by*px) * sign(-t._pr2)
    if nc < -_abs_err_intriangle_zero
        return -3
    elseif nc < _abs_err_intriangle_zero
        return _sz_intriangle(t, p, px, py)
    end

    const l = nb+nc + t._pr2*sign(-t._pr2)
    if l < -_abs_err_intriangle
        return 1
    elseif l > _abs_err_intriangle
        return -1
    end
    _sz_intriangle(t, p, px, py)
end

# fine filtered intriangle
function _sz_intriangle(t::TetrahedronTypes, p::AbstractPoint3D)
    const px = getx(p) - getx(geta(t)); const py = gety(p) - gety(geta(t)); const pz = getz(p) - getz(geta(t))
    const sz = abs(px) + abs(py) + abs(pz) +
                    abs(t._bx) + abs(t._by) + abs(t._bz) +
                    abs(t._cx) + abs(t._cy) + abs(t._cz) +
                    abs(t._dx) + abs(t._dy) + abs(t._dz)

    const nb = (t._cx*t._dy*pz-t._cx*t._dz*py-t._cy*t._dx*pz+t._cy*t._dz*px+t._cz*t._dx*py-t._cz*t._dy*px) * sign(-t._pr2)
    if nb < -_abs_err_intetra_zero*sz
        return -2
    elseif nb < _abs_err_intetra_zero*sz
        # we need an exact calculation
        return _exact_intriangle(t, p)
    end

    const nc = (-t._bx*t._dy*pz+t._bx*t._dz*py+t._by*t._dx*pz-t._by*t._dz*px-t._bz*t._dx*py+t._bz*t._dy*px) * sign(-t._pr2)
    if nc < -_abs_err_intetra_zero*sz
        return -3
    elseif nc < _abs_err_intetra_zero*sz
        # we need an exact calculation
        return _exact_intriangle(t, p)
    end

    const nd = (t._bx*t._cy*pz-t._bx*t._cz*py-t._by*t._cx*pz+t._by*t._cz*px+t._bz*t._cx*py-t._bz*t._cy*px) * sign(-t._pr2)
    if nd < -_abs_err_intetra_zero*sz
        return -4
    elseif nd < _abs_err_intetra_zero*sz
        # we need an exact calculation
        return _exact_intriangle(t, p)
    end

    const l = nb+nc+nd + t._pr2*sign(-t._pr2)
    if l < -_abs_err_intetra*sz
        return 1
    elseif l > _abs_err_intetra*sz
        return -1
    end
    # we need an exact calculation
    _exact_intriangle(t, p)
end

# gross filtered intriangle
function intriangle(t::TetrahedronTypes, p::AbstractPoint3D)
    const px = getx(p) - getx(geta(t)); const py = gety(p) - gety(geta(t)); const pz = getz(p) - getz(geta(t))

    const nb = (t._cx*t._dy*pz-t._cx*t._dz*py-t._cy*t._dx*pz+t._cy*t._dz*px+t._cz*t._dx*py-t._cz*t._dy*px) * sign(-t._pr2)
    if nb < -_abs_err_intetra_zero
        return -2
    elseif nb < _abs_err_intetra_zero
        return _sz_intriangle(t, p)
    end

    const nc = (-t._bx*t._dy*pz+t._bx*t._dz*py+t._by*t._dx*pz-t._by*t._dz*px-t._bz*t._dx*py+t._bz*t._dy*px) * sign(-t._pr2)
    if nc < -_abs_err_intetra_zero
        return -3
    elseif nc < _abs_err_intetra_zero
        return _sz_intriangle(t, p)
    end

    const nd = (t._bx*t._cy*pz-t._bx*t._cz*py-t._by*t._cx*pz+t._by*t._cz*px+t._bz*t._cx*py-t._bz*t._cy*px) * sign(-t._pr2)
    if nd < -_abs_err_intetra_zero
        return -4
    elseif nd < _abs_err_intetra_zero
        return _sz_intriangle(t, p)
    end

    const l = nb+nc+nd + t._pr2*sign(-t._pr2)
    if l < -_abs_err_intetra
        return 1
    elseif l > _abs_err_intetra
        return -1
    end
    _sz_intriangle(t, p)
end

# helper  methods to use the filtered (fast, exact) intriangle directly with raw coordinates
intriangle(ax::Float64, ay::Float64, bx::Float64, by::Float64, cx::Float64, cy::Float64, px::Float64, py::Float64) =
    intriangle(Triangle(ax, ay, bx, by, cx, cy), Point2D(px, py))
intriangle(ax::Float64, ay::Float64, az::Float64, bx::Float64, by::Float64, bz::Float64, cx::Float64, cy::Float64, cz::Float64, dx::Float64, dy::Float64, dz::Float64, px::Float64, py::Float64, pz::Float64) =
    intriangle(Tetrahedron(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz), Point3D(px, py, pz))


###################################################################################
#
#    Hilbert stuff
#
#

# number of bits to use per dimension in calculating the peano-key
const peano_2D_bits = 31
const peano_3D_bits = 21

 # implementing 2D scale dependednt Peano-Hilbert indexing

if VERSION < v"0.4-"
    _extract_peano_bin_num(nbins::Int64, n::Float64) = itrunc(Integer, (n-1)*nbins )
else
    _extract_peano_bin_num(nbins::Int64, n::Float64) = trunc(Integer, (n-1)*nbins )
end

# calculate peano key for given point
function peanokey(p::AbstractPoint2D, bits::Int64=peano_2D_bits)
    const n = 1 << bits
    s = n >> 1; d = 0
    x = _extract_peano_bin_num(n, getx(p))
    y = _extract_peano_bin_num(n, gety(p))
    while true
        rx = (x & s) > 0
        ry = (y & s) > 0
        d += s * s * ((3 * rx) $ ry)
        s = s >> 1
        (s == 0) && break
        if ry == 0
            if rx == 1
                x = n - 1 - x;
                y = n - 1 - y;
            end
            x, y = y, x
        end
    end
    d
end

# Inverse calculation. I.e. calculate the point that given given peano key
function Point2D(peanokey::Int64, bits::Int64=peano_2D_bits)
    const n = 1 << bits
    x = 0; y = 0; s=1
    while true

        rx = 1 & (peanokey >> 1)
        ry = 1 & (peanokey $ rx)

        if ry == 0
            if rx == 1
                x = s - 1 - x;
                y = s - 1 - y;
            end
            x, y = y, x
        end

        x += s * rx
        y += s * ry

        s = s << 1
        (s >= n) && break

        peanokey = peanokey >> 2
    end
    Point2D(1+x/n, 1+y/n)
end

# implementing 3D scaleful Peano-Hilbert indexing

const quadrants_arr = [
  0, 7, 1, 6, 3, 4, 2, 5,
  7, 4, 6, 5, 0, 3, 1, 2,
  4, 3, 5, 2, 7, 0, 6, 1,
  3, 0, 2, 1, 4, 7, 5, 6,
  1, 0, 6, 7, 2, 3, 5, 4,
  0, 3, 7, 4, 1, 2, 6, 5,
  3, 2, 4, 5, 0, 1, 7, 6,
  2, 1, 5, 6, 3, 0, 4, 7,
  6, 1, 7, 0, 5, 2, 4, 3,
  1, 2, 0, 3, 6, 5, 7, 4,
  2, 5, 3, 4, 1, 6, 0, 7,
  5, 6, 4, 7, 2, 1, 3, 0,
  7, 6, 0, 1, 4, 5, 3, 2,
  6, 5, 1, 2, 7, 4, 0, 3,
  5, 4, 2, 3, 6, 7, 1, 0,
  4, 7, 3, 0, 5, 6, 2, 1,
  6, 7, 5, 4, 1, 0, 2, 3,
  7, 0, 4, 3, 6, 1, 5, 2,
  0, 1, 3, 2, 7, 6, 4, 5,
  1, 6, 2, 5, 0, 7, 3, 4,
  2, 3, 1, 0, 5, 4, 6, 7,
  3, 4, 0, 7, 2, 5, 1, 6,
  4, 5, 7, 6, 3, 2, 0, 1,
  5, 2, 6, 1, 4, 3, 7, 0]
quadrants(a::Int64, b::Int64, c::Int64, d::Int64) = (@inbounds const x = quadrants_arr[1+a<<3+b<<2+c<<1+d]; x)
rotxmap_table = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22]
rotymap_table = [1, 2, 3, 0, 16, 17, 18, 19, 11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7]
rotx_table = [3, 0, 0, 2, 2, 0, 0, 1]
roty_table = [0, 1, 1, 2, 2, 3, 3, 0]
sense_table = [-1, -1, -1, +1, +1, -1, -1, -1]

function peanokey(p::AbstractPoint3D, bits::Int64=peano_3D_bits)
    const n = 1 << bits
    const x = _extract_peano_bin_num(n, getx(p))
    const y = _extract_peano_bin_num(n, gety(p))
    const z = _extract_peano_bin_num(n, getz(p))
    mask = 1 << (bits - 1)
    key = 0
    rotation = 0
    sense = 1
    for i in 1:bits
        const bitx = (x & mask > 0) ? 1 : 0
        const bity = (y & mask > 0) ? 1 : 0
        const bitz = (z & mask > 0) ? 1 : 0

        const quad = quadrants(rotation, bitx, bity, bitz)

        key <<= 3
        key += sense == 1 ? quad : 7-quad

        @inbounds rotx = rotx_table[quad+1]
        @inbounds roty = roty_table[quad+1]
        @inbounds sense *= sense_table[quad+1]

        while rotx > 0
            @inbounds rotation = rotxmap_table[rotation+1]
            rotx -= 1
        end

        while(roty > 0)
            @inbounds rotation = rotymap_table[rotation+1]
            roty -= 1
        end
        mask >>= 1
    end

    key
end

const quadrants_inverse_x = Array(Int64, (24,8))
const quadrants_inverse_y = Array(Int64, (24,8))
const quadrants_inverse_z = Array(Int64, (24,8))

function _init_inv_peano_3d()
    for rotation in 0:23
        for bitx in 0:1
            for bity in 0:1
                for bitz in 0:1
                    quad = quadrants(rotation, bitx, bity, bitz)
                    quadrants_inverse_x[rotation+1, quad+1] = bitx;
                    quadrants_inverse_y[rotation+1, quad+1] = bity;
                    quadrants_inverse_z[rotation+1, quad+1] = bitz;
                end
            end
        end
    end
end

_init_inv_peano_3d()

# inverse transformation, give a key get a point
function Point3D(key::Int64, bits::Int64=peano_3D_bits)

    const n = 1 << bits

    shift = 3*(bits - 1)
    mask = 7 << shift

    rotation = 0
    sense = 1

    x = 0
    y = 0
    z = 0

    for i in 1:bits
        keypart = (key & mask) >> shift

        quad = sense == 1 ? keypart : 7 - keypart

        x = (x << 1) + quadrants_inverse_x[rotation+1, quad+1]
        y = (y << 1) + quadrants_inverse_y[rotation+1, quad+1]
        z = (z << 1) + quadrants_inverse_z[rotation+1, quad+1]

        rotx = rotx_table[quad+1]
        roty = roty_table[quad+1]
        sense *= sense_table[quad+1]

        while rotx > 0
            rotation = rotxmap_table[rotation+1]
            rotx -= 1
        end

        while roty > 0
            rotation = rotymap_table[rotation+1]
            roty -= 1
        end

        mask >>= 3
        shift -= 3
    end
    Point3D(1+x/n, 1+y/n, 1+z/n)
end

# implementing scale-free Hilbert ordering. Real all about it here:
# http://doc.cgal.org/latest/Spatial_sorting/index.html

abstract AbstractCoordinate
type CoordinateX <: AbstractCoordinate end
type CoordinateY <: AbstractCoordinate end
type CoordinateZ <: AbstractCoordinate end
const coordinatex = CoordinateX()
const coordinatey = CoordinateY()
const coordinatez = CoordinateZ()
next2d(::CoordinateX) = coordinatey
next2d(::CoordinateY) = coordinatex
next3d(::CoordinateX) = coordinatey
next3d(::CoordinateY) = coordinatez
next3d(::CoordinateZ) = coordinatex
nextnext3d(::CoordinateX) = coordinatez
nextnext3d(::CoordinateY) = coordinatex
nextnext3d(::CoordinateZ) = coordinatey

abstract AbstractDirection
type Forward <: AbstractDirection end
type Backward <: AbstractDirection end
const forward = Forward()
const backward = Backward()
@compat Base.:!(::Forward) = backward
@compat Base.:!(::Backward) = forward

compare(::Forward, ::CoordinateX, p1::AbstractPoint, p2::AbstractPoint) = getx(p1) < getx(p2)
compare(::Backward, ::CoordinateX, p1::AbstractPoint, p2::AbstractPoint) = getx(p1) > getx(p2)
compare(::Forward, ::CoordinateY, p1::AbstractPoint, p2::AbstractPoint) = gety(p1) < gety(p2)
compare(::Backward, ::CoordinateY, p1::AbstractPoint, p2::AbstractPoint) = gety(p1) > gety(p2)
compare(::Forward, ::CoordinateZ, p1::AbstractPoint, p2::AbstractPoint) = getz(p1) < getz(p2)
compare(::Backward, ::CoordinateZ, p1::AbstractPoint, p2::AbstractPoint) = getz(p1) > getz(p2)

function select!{T<:AbstractPoint}(direction::AbstractDirection, coordinate::AbstractCoordinate, v::Array{T,1}, k::Int, lo::Int, hi::Int)
    lo <= k <= hi || error("select index $k is out of range $lo:$hi")
    @inbounds while lo < hi
        if hi-lo == 1
            if compare(direction, coordinate, v[hi], v[lo])
                v[lo], v[hi] = v[hi], v[lo]
            end
            return v[k]
        end
        pivot = v[(lo+hi)>>>1]
        i, j = lo, hi
        while true
            while compare(direction, coordinate, v[i], pivot); i += 1; end
            while compare(direction, coordinate, pivot, v[j]); j -= 1; end
            i <= j || break
            v[i], v[j] = v[j], v[i]
            i += 1; j -= 1
        end
        if k <= j
            hi = j
        elseif i <= k
            lo = i
        else
            return pivot
        end
    end
    return v[lo]
end

function hilbertsort!{T<:AbstractPoint2D}(directionx::AbstractDirection, directiony::AbstractDirection, coordinate::AbstractCoordinate, a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64=4)
    hi-lo <= lim && return a

    const i2 = (lo+hi)>>>1
    const i1 = (lo+i2)>>>1
    const i3 = (i2+hi)>>>1

    select!(directionx, coordinate, a, i2, lo, hi)
    select!(directiony, next2d(coordinate), a, i1, lo, i2)
    select!(!directiony, next2d(coordinate), a, i3, i2, hi)

    hilbertsort!(directiony, directionx, next2d(coordinate), a, lo, i1, lim)
    hilbertsort!(directionx, directiony, coordinate, a, i1, i2, lim)
    hilbertsort!(directionx, directiony, coordinate, a, i2, i3, lim)
    hilbertsort!(!directiony, !directionx, next2d(coordinate), a, i3, hi, lim)

    return a
end

function hilbertsort!{T<:AbstractPoint3D}(directionx::AbstractDirection, directiony::AbstractDirection, directionz::AbstractDirection, coordinate::AbstractCoordinate, a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64=8)
    hi-lo <= lim && return a

    const i4 = (lo+hi)>>>1
    const i2 = (lo+i4)>>>1
    const i1 = (lo+i2)>>>1
    const i3 = (i2+i4)>>>1
    const i6 = (i4+hi)>>>1
    const i5 = (i4+i6)>>>1
    const i7 = (i6+hi)>>>1

    select!(directionx, coordinate, a, i4, lo, hi)
    select!(directiony, next3d(coordinate), a, i2, lo, i4)
    select!(directionz, nextnext3d(coordinate), a, i1, lo, i2)
    select!(!directionz, nextnext3d(coordinate), a, i3, i2, i4)
    select!(!directiony, next3d(coordinate), a, i6, i4, hi)
    select!(directionz, nextnext3d(coordinate), a, i5, i4, i6)
    select!(!directionz, nextnext3d(coordinate), a, i7, i6, hi)

    hilbertsort!( directionz,  directionx,  directiony, nextnext3d(coordinate), a, lo, i1, lim)
    hilbertsort!( directiony,  directionz,  directionx, next3d(coordinate),     a, i1, i2, lim)
    hilbertsort!( directiony,  directionz,  directionx, next3d(coordinate),     a, i2, i3, lim)
    hilbertsort!( directionx, !directiony, !directionz, coordinate,             a, i3, i4, lim)
    hilbertsort!( directionx, !directiony, !directionz, coordinate,             a, i4, i5, lim)
    hilbertsort!(!directiony,  directionz, !directionx, next3d(coordinate),     a, i5, i6, lim)
    hilbertsort!(!directiony,  directionz, !directionx, next3d(coordinate),     a, i6, i7, lim)
    hilbertsort!(!directionz, !directionx,  directiony, nextnext3d(coordinate), a, i7, hi, lim)

    return a
end

hilbertsort!{T<:AbstractPoint2D}(a::Array{T,1}) = hilbertsort!(backward, backward, coordinatey, a, 1, length(a))
hilbertsort!{T<:AbstractPoint2D}(a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64) = hilbertsort!(backward, backward, coordinatey, a, lo, hi, lim)
hilbertsort!{T<:AbstractPoint3D}(a::Array{T,1}) = hilbertsort!(backward, backward, backward, coordinatez, a, 1, length(a))
hilbertsort!{T<:AbstractPoint3D}(a::Array{T,1}, lo::Int64, hi::Int64, lim::Int64) = hilbertsort!(backward, backward, backward, coordinatey, a, lo, hi, lim)

# multi-scale sort. Read all about it here:
# http://doc.cgal.org/latest/Spatial_sorting/classCGAL_1_1Multiscale__sort.html
function _mssort!{T<:AbstractPoint}(a::Array{T,1}, lim_ms::Int64, lim_hl::Int64, rat::Float64)
    hi = length(a)
    lo = 1
    while true
        lo = hi - round(Int, (1-rat)*hi)
        hi-lo <= lim_ms && return a
        hilbertsort!(a, lo, hi, lim_hl)
        hi = lo-1
    end
    return a
end

# Utility methods, setting some different defaults for 2D and 3D. These are exported
mssort!{T<:AbstractPoint2D}(a::Array{T,1}; lim_ms::Int64=16, lim_hl::Int64=4, rat::Float64=0.25) =
    _mssort!(a, lim_ms, lim_hl, rat)
mssort!{T<:AbstractPoint3D}(a::Array{T,1}; lim_ms::Int64=64, lim_hl::Int64=8, rat::Float64=0.125) =
    _mssort!(a, lim_ms, lim_hl, rat)

end # module GeometicalPredicates
