using GeometricalPredicates
using Test

@testset "2D orientation" begin
    @testset "Lines" begin
        a = Point(1.1, 1.1)
        b = Point(1.5, 1.5)

        l = Line(a, b)
        @test orientation(l, Point(1.4, 1.6)) == 1
        @test orientation(l, Point(1.4,1.05)) == -1
        @test orientation(l, Point(1.4,1.4)) == 0

        l = Line(b, a)
        @test orientation(l, Point(1.4, 1.6)) == -1
        @test orientation(l, Point(1.4,1.05)) == 1
        @test orientation(l, Point(1.4,1.4)) == 0

        @test abs(length2(l)-0.4*0.4*2) < 1e-5
    end

    @testset "Triangles" begin
        ax = 1.1; ay = 1.1
        bx = 1.2; by = 1.2
        cx = 1.3; cy = 1.3
        @test orientation(ax,ay,bx,by,cx,cy) == 0
        @test orientation(bx,by,ax,ay,cx,cy) == 0
        @test orientation(bx,by,cx,cy,ax,ay) == 0

        ax = 1.1; ay = 1.1
        bx = 1.2; by = 1.199999999999
        cx = 1.3; cy = 1.3
        @test orientation(ax,ay,bx,by,cx,cy) == 1
        @test orientation(bx,by,ax,ay,cx,cy) == -1
        @test orientation(bx,by,cx,cy,ax,ay) == 1

        ax = 1.1; ay = 1.1
        bx = 1.2; by = 1.199999999999999
        cx = 1.3; cy = 1.3
        @test orientation(ax,ay,bx,by,cx,cy) == 1
        @test orientation(bx,by,ax,ay,cx,cy) == -1
        @test orientation(bx,by,cx,cy,ax,ay) == 1

        ax = 1.1; ay = 1.1
        bx = 1.2; by = 1.21
        cx = 1.3; cy = 1.3
        @test orientation(ax,ay,bx,by,cx,cy) == -1
        @test orientation(bx,by,ax,ay,cx,cy) == 1
        @test orientation(bx,by,cx,cy,ax,ay) == -1

        ax = 1.1; ay = 1.1
        bx = 1.2; by = 1.20000000000001
        cx = 1.3; cy = 1.3
        @test orientation(ax,ay,bx,by,cx,cy) == -1
        @test orientation(bx,by,ax,ay,cx,cy) == 1
        @test orientation(bx,by,cx,cy,ax,ay) == -1

        p1 = Point2D(1.0, 1.0)
        p2 = Point2D(1.9, 1.5)
        p3 = Point2D(1.45, 1.25)
        tr = Triangle(p1, p2, p3)
        @test orientation(tr) == 0

        p3 = Point2D(getx(p3), 1.3)
        tr = Triangle(p1, p2, p3)
        @test orientation(tr) == 1

        p3 = Point2D(getx(p3), 1.2)
        tr = Triangle(p1, p2, p3)
        @test orientation(tr) == -1

        tr = Triangle(p1, p2, p3, positivelyoriented)
        @test orientation(tr) == 1

        tr = Triangle(p1, p3, p2, negativelyoriented)
        @test orientation(tr) == -1
    end
end

@testset "2D incircle" begin
    ax = 1.1; ay = 1.1
    bx = 1.2; by = 1.2
    cx = 1.3; cy = 1.3
    dx = 1.4; dy = 1.7
    @test incircle(ax,ay,bx,by,cx,cy,dx,dy) == 2

    ax = 1.1; ay = 1.1
    bx = 1.3; by = 1.1
    cx = 1.3; cy = 1.3
    dx = 1.1; dy = 1.3
    @test incircle(ax,ay,bx,by,cx,cy,dx,dy) == 0
    @test incircle(bx,by,ax,ay,cx,cy,dx,dy) == 0
    @test incircle(bx,by,cx,cy,ax,ay,dx,dy) == 0

    ax = 1.1; ay = 1.1
    bx = 1.3; by = 1.1
    cx = 1.3; cy = 1.3
    dx = 1.1; dy = 1.3000001
    @test incircle(ax,ay,bx,by,cx,cy,dx,dy) == -1
    @test incircle(bx,by,ax,ay,cx,cy,dx,dy) == -1
    @test incircle(bx,by,cx,cy,ax,ay,dx,dy) == -1

    ax = 1.1; ay = 1.1
    bx = 1.3; by = 1.1
    cx = 1.3; cy = 1.3
    dx = 1.1; dy = 1.30000000000001
    @test incircle(ax,ay,bx,by,cx,cy,dx,dy) == -1
    @test incircle(bx,by,ax,ay,cx,cy,dx,dy) == -1
    @test incircle(bx,by,cx,cy,ax,ay,dx,dy) == -1

    ax = 1.1; ay = 1.1
    bx = 1.3; by = 1.1
    cx = 1.3; cy = 1.3
    dx = 1.1; dy = 1.29999
    @test incircle(ax,ay,bx,by,cx,cy,dx,dy) == 1
    @test incircle(bx,by,ax,ay,cx,cy,dx,dy) == 1
    @test incircle(bx,by,cx,cy,ax,ay,dx,dy) == 1

    ax = 1.1; ay = 1.1
    bx = 1.3; by = 1.1
    cx = 1.3; cy = 1.3
    dx = 1.1; dy = 1.29999999999999
    @test incircle(bx,by,ax,ay,cx,cy,dx,dy) == 1
    @test incircle(bx,by,ax,ay,cx,cy,dx,dy) == 1
    @test incircle(bx,by,cx,cy,ax,ay,dx,dy) == 1

    p1 = Point2D(1.0, 1.0)
    p2 = Point2D(1.0625, 1.0)
    p3 = Point2D(1.0625, 1.0625)
    tr = Triangle(p1, p2, p3)

    p4 = Point2D(1.03, 1.003)
    p5 = Point2D(1.99, 1.99)
    p6 = Point2D(1.0, 1.0625)
    @test incircle(tr, p4) == 1
    @test incircle(tr, p5) == -1
    @test incircle(tr, p6) == 0
    @test incircle(tr, p2) == 0
    @test incircle(tr, p3) == 0
    @test incircle(tr, p1) == 0
end

@testset "3D orientation" begin
    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.2; by = 1.2; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.4; dy = 1.4; dz=1.1
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == 0
    @test orientation(dx,dy,dz,bx,by,bz,cx,cy,cz,ax,ay,az) == 0
    @test orientation(dx,dy,dz,cx,cy,cz,bx,by,bz,ax,ay,az) == 0

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.2; by = 1.2; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.4; dy = 1.4; dz=1.100001
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == 0
    @test orientation(dx,dy,dz,bx,by,bz,cx,cy,cz,ax,ay,az) == 0
    @test orientation(dx,dy,dz,cx,cy,cz,bx,by,bz,ax,ay,az) == 0

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.2; by = 1.2; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.4; dy = 1.4; dz=1.100000000001
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == 0
    @test orientation(dx,dy,dz,bx,by,bz,cx,cy,cz,ax,ay,az) == 0
    @test orientation(dx,dy,dz,cx,cy,cz,bx,by,bz,ax,ay,az) == 0

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.2; by = 1.2; bz=1.1
    cx = 1.3; cy = 1.15; cz=1.1
    dx = 1.4; dy = 1.17; dz=1.1
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == 0
    @test orientation(dx,dy,dz,bx,by,bz,cx,cy,cz,ax,ay,az) == 0
    @test orientation(dx,dy,dz,cx,cy,cz,bx,by,bz,ax,ay,az) == 0

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.2; dy = 1.2; dz=1.5
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == 1
    @test orientation(dx,dy,dz,bx,by,bz,cx,cy,cz,ax,ay,az) == -1
    @test orientation(dx,dy,dz,cx,cy,cz,bx,by,bz,ax,ay,az) == 1

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.2; dy = 1.2; dz=1.10000000000001
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == 1
    @test orientation(bx,by,bz,ax,ay,az,cx,cy,cz,dx,dy,dz) == -1

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.2; dy = 1.2; dz=1.09
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == -1
    @test orientation(bx,by,bz,ax,ay,az,cx,cy,cz,dx,dy,dz) == 1
    @test orientation(ax,ay,az,bx,by,bz,dx,dy,dz,cx,cy,cz) == 1

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.2; dy = 1.2; dz=1.099999999999999
    @test orientation(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz) == -1
    @test orientation(bx,by,bz,ax,ay,az,cx,cy,cz,dx,dy,dz) == 1
    @test orientation(ax,ay,az,bx,by,bz,dx,dy,dz,cx,cy,cz) == 1

    p1 = Point3D(ax, ay, az)
    p2 = Point3D(bx, by, bz)
    p3 = Point3D(cx, cy, cz)
    p4 = Point3D(dx, dy, dz)

    t = Tetrahedron(p1, p2, p3, p4)
    @test orientation(t) == -1

    t = Tetrahedron(p1, p2, p3, p4, positivelyoriented)
    @test orientation(t) == 1

    t = Tetrahedron(p1, p2, p4, p3, negativelyoriented)
    @test orientation(t) == -1
end

@testset "3D incircle" begin
    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.4; dy = 1.4; dz=1.2
    ex = 1.1; ey = 1.1; ez=1.1
    @test incircle(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez) == 0
    @test incircle(bx,by,bz,ax,ay,az,cx,cy,cz,dx,dy,dz,ex,ey,ez) == 0
    @test incircle(ax,ay,az,bx,by,bz,dx,dy,dz,cx,cy,cz,ex,ey,ez) == 0

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.2; dy = 1.2; dz=1.5
    ex = 1.2; ey = 1.15; ez=1.13
    @test incircle(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez) == 1
    @test incircle(bx,by,bz,ax,ay,az,cx,cy,cz,dx,dy,dz,ex,ey,ez) == 1
    @test incircle(ax,ay,az,bx,by,bz,dx,dy,dz,cx,cy,cz,ex,ey,ez) == 1

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.2; dy = 1.2; dz=1.5
    ex = 1.2; ey = 1.15; ez=1.1000000000001
    @test incircle(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez) == 1
    @test incircle(bx,by,bz,ax,ay,az,cx,cy,cz,dx,dy,dz,ex,ey,ez) == 1
    @test incircle(ax,ay,az,bx,by,bz,dx,dy,dz,cx,cy,cz,ex,ey,ez) == 1

    ax = 1.1; ay = 1.1; az=1.1
    bx = 1.3; by = 1.1; bz=1.1
    cx = 1.3; cy = 1.3; cz=1.1
    dx = 1.2; dy = 1.2; dz=1.5
    ex = 1.9; ey = 1.15; ez=1.01
    @test incircle(ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,ex,ey,ez) == -1
    @test incircle(bx,by,bz,ax,ay,az,cx,cy,cz,dx,dy,dz,ex,ey,ez) == -1
    @test incircle(ax,ay,az,bx,by,bz,dx,dy,dz,cx,cy,cz,ex,ey,ez) == -1
end

@testset "intriangle (2D)!" begin
    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.11; dy = 1.11
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == 1
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) == 1
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) == 1

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.01; dy = 1.01
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) < 0
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) < 0
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) < 0

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.31; dy = 1.01
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) < 0
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) < 0
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) < 0

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.1; dy = 1.2
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == 3
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) == 2
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) == 2

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.2; dy = 1.1
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == 4
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) == 4
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) == 3

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.25; dy = 1.25
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == 2
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) == 3
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) == 4

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.25; dy = 1.250000000000001
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) < 0
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) < 0
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) < 0

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.25; dy = 1.249999999999999
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == 1
    @test intriangle(bx,by,ax,ay,cx,cy,dx,dy) == 1
    @test intriangle(bx,by,cx,cy,ax,ay,dx,dy) == 1

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.7; dy = 1.1
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -1
    dx = 1.1; dy = 1.7
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -1
    dx = 1.01; dy = 1.1
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -2
    dx = 1.1; dy = 1.01
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -3

    ax = 1.1; ay = 1.1
    bx = 1.4; by = 1.1
    cx = 1.1; cy = 1.4
    dx = 1.55; dy = 1.55
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -1
    dx = 1.09; dy = 1.11
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -2
    dx = 1.11; dy = 1.09
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -3

    ax = min_coord; ay = min_coord
    bx = 1.1; by = 1.1
    cx = max_coord; cy = min_coord
    dx = 1.2; dy = 1.2
    @test intriangle(ax,ay,bx,by,cx,cy,dx,dy) == -1
end

@testset "intriangle (3D) !" begin
    ax = 1.1; ay = 1.1; az = 1.1
    bx = 1.4; by = 1.1; bz = 1.1
    cx = 1.1; cy = 1.4; cz = 1.1
    dx = 1.1; dy = 1.1; dz = 1.4
    px = 1.11; py = 1.11; pz = 1.11
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == 1
    @test intriangle(bx, by, bz, ax, ay, az, cx, cy, cz, dx, dy, dz, px, py, pz) == 1
    @test intriangle(ax, ay, az, cx, cy, cz, bx, by, bz, dx, dy, dz, px, py, pz) == 1

    ax = 1.1; ay = 1.1; az = 1.1
    bx = 1.4; by = 1.1; bz = 1.1
    cx = 1.1; cy = 1.4; cz = 1.1
    dx = 1.1; dy = 1.1; dz = 1.4
    px = 1.02; py = 1.2; pz = 1.2
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) < 0
    @test intriangle(bx, by, bz, ax, ay, az, cx, cy, cz, dx, dy, dz, px, py, pz) < 0
    @test intriangle(ax, ay, az, cx, cy, cz, bx, by, bz, dx, dy, dz, px, py, pz) < 0

    ax = 1.1; ay = 1.1; az = 1.1
    bx = 1.4; by = 1.1; bz = 1.1
    cx = 1.1; cy = 1.4; cz = 1.1
    dx = 1.1; dy = 1.1; dz = 1.4
    px = 1.09999999999999; py = 1.2; pz = 1.2
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) < 0
    @test intriangle(bx, by, bz, ax, ay, az, cx, cy, cz, dx, dy, dz, px, py, pz) < 0
    @test intriangle(ax, ay, az, cx, cy, cz, bx, by, bz, dx, dy, dz, px, py, pz) < 0

    ax = 1.1; ay = 1.1; az = 1.1
    bx = 1.4; by = 1.1; bz = 1.1
    cx = 1.1; cy = 1.4; cz = 1.1
    dx = 1.1; dy = 1.1; dz = 1.4
    px = 1.10000000000001; py = 1.2; pz = 1.2
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == 1
    @test intriangle(bx, by, bz, ax, ay, az, cx, cy, cz, dx, dy, dz, px, py, pz) == 1
    @test intriangle(ax, ay, az, cx, cy, cz, bx, by, bz, dx, dy, dz, px, py, pz) == 1

    ax = 1.1; ay = 1.1; az = 1.1
    bx = 1.4; by = 1.1; bz = 1.1
    cx = 1.1; cy = 1.4; cz = 1.1
    dx = 1.1; dy = 1.1; dz = 1.4
    px = 1.11; py = 1.11; pz = 1.1
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == 5
    @test intriangle(bx, by, bz, ax, ay, az, cx, cy, cz, dx, dy, dz, px, py, pz) == 5
    @test intriangle(ax, ay, az, cx, cy, cz, bx, by, bz, dx, dy, dz, px, py, pz) == 5
    @test intriangle(dx, dy, dz, bx, by, bz, cx, cy, cz, ax, ay, az, px, py, pz) == 2
    @test intriangle(ax, ay, az, dx, dy, dz, cx, cy, cz, bx, by, bz, px, py, pz) == 3
    @test intriangle(ax, ay, az, bx, by, bz, dx, dy, dz, cx, cy, cz, px, py, pz) == 4

    ax = 1.1; ay = 1.1; az = 1.1
    bx = 1.4; by = 1.1; bz = 1.1
    cx = 1.1; cy = 1.4; cz = 1.1
    dx = 1.1; dy = 1.1; dz = 1.4
    px = 1.55; py = 1.55; pz = 1.55
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.09; py = 1.11; pz = 1.11
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -2
    px = 1.11; py = 1.09; pz = 1.11
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -3
    px = 1.11; py = 1.11; pz = 1.09
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -4

    ax = 1.1; ay = 1.1; az = 1.1
    bx = 1.4; by = 1.1; bz = 1.1
    cx = 1.1; cy = 1.4; cz = 1.1
    dx = 1.1; dy = 1.1; dz = 1.4
    px = 1.7; py = 1.1; pz = 1.1
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.7; py = 1.1; pz = 1.1
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.1; py = 1.7; pz = 1.1
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.1; py = 1.1; pz = 1.7
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.9; py = 1.9; pz = 1.1
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.9; py = 1.1; pz = 1.9
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.1; py = 1.9; pz = 1.9
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -1
    px = 1.01; py = 1.1; pz = 1.1
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -2
    px = 1.1; py = 1.01; pz = 1.1
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -3
    px = 1.1; py = 1.1; pz = 1.01
    @test intriangle(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz, px, py, pz) == -4
end

@testset "Scale dependednt Peano-Hilbert stuff" begin
    @testset "2D" begin
        @test peanokey(Point2D(1.01,1.01), 2) == 0
        @test peanokey(Point2D(1.9901,1.01), 2) == 4*4-1
        @test Point2D(15, 2) == Point2D(1.75, 1.0)
        @test Point2D(0, 2) == Point2D(1.0, 1.0)
        for x in range(1.0,stop=1.999999,length=100), y in range(1.0,stop=1.999999,length=100)
            p = Point2D(x, y)
            d = peanokey(p)
            pp= Point2D(d)
            @test abs(getx(p) - getx(pp)) < 1e-7
            @test abs(gety(p) - gety(pp)) < 1e-7
        end
    end

    @testset "3D" begin
        @test peanokey(Point3D(1.01,1.01,1.01), 2) == 0
        @test peanokey(Point3D(1.01,1.01,1.9901), 2) == 4*4*4-1
        @test Point3D(15, 2) == Point3D(1.25,1.5,1.0)
        @test Point3D(0, 2) == Point3D(1.0, 1.0, 1.0)
        for x in range(1.0,stop=1.99999,length=30), y in range(1.0,stop=1.99999,length=30), z in range(1.0,stop=1.999999,length=30)
            p = Point3D(x, y, z)
            d = peanokey(p)
            pp= Point3D(d)
            @test abs(getx(p) - getx(pp)) < 1e-5
            @test abs(gety(p) - gety(pp)) < 1e-5
            @test abs(getz(p) - getz(pp)) < 1e-5
        end
    end
end

@testset "Scale independent Peano-Hilbert stuff" begin
    a2 = Point2D[]
    N=1000
    for i in 1:N
        push!(a2, Point2D(1.0+rand(), 1.0+rand()))
    end
    # TODO: not much of a test, I need to betteer think how to test this
    # At least this runs the code...
    @test length(mssort!(a2)) == N

    a3 = Point3D[]
    for i in 1:N
        push!(a3, Point3D(1.0+rand(), 1.0+rand(), 1.0+rand()))
    end
    # TODO: not much of a test, I need to betteer think how to test this
    # At least this runs the code...
    @test length(mssort!(a3)) == N
end

@testset "Qttys functions" begin
    @testset "2D" begin
        a = Point2D(1.0, 1.0)
        b = Point2D(1.0, 1.5)
        c = Point2D(1.5, 1.0)
        t = Triangle(a, b, c)
        @test abs(area(t)-1/8) < 1e-7
        c = centroid(t)
        @test abs(getx(c)-1-1/6) < 1e-7
        @test abs(gety(c)-1-1/6) < 1e-7
        c = circumcenter(t)
        @test abs(getx(c)-1.25) < 1e-7
        @test abs(gety(c)-1.25) < 1e-7
        @test abs(circumradius2(t)-1/8) < 1e-7
    end

    @testset "3D" begin
        a = Point3D(1.0, 1.0, 1.0)
        b = Point3D(1.0, 1.0, 1.5)
        c = Point3D(1.0, 1.5, 1.0)
        d = Point3D(1.5, 1.0, 1.0)
        t = Tetrahedron(a, b, c, d)
        @test abs(volume(t)-1/16) < 1e-7
        c = centroid(t)
        @test abs(getx(c)-1.125) < 1e-7
        @test abs(gety(c)-1.125) < 1e-7
        @test abs(getz(c)-1.125) < 1e-7
        c = circumcenter(t)
        @test abs(getx(c)-1.25) < 1e-7
        @test abs(gety(c)-1.25) < 1e-7
        @test abs(getz(c)-1.25) < 1e-7
        @test abs(circumradius2(t)-0.25^2*3) < 1e-7
    end
end

@testset "2D Polygons" begin
    ll = Point(1.0,1.0)
    lr = Point(1.2,1.0)
    ur = Point(1.2,1.2)
    ul = Point(1.0,1.2)
    poly = Polygon(ll, lr, ur, ul)
    @test isempty( setdiff(getpoints(poly), [ll,lr,ur,ul]) )
    @test isempty( setdiff(getlines(poly), [Line(ll,lr), Line(lr,ur), Line(ur,ul), Line(ul,ll)]) )
    @test inpolygon(poly, Point(1.1,1.1)) == true
    @test inpolygon(poly, Point(1.1,1.3)) == false
end
# that's it for today
