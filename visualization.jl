using Makie, CDDLib, Polyhedra

v = convexhull([0, 0, 0]) + conichull([1, 0, 0], [0, 1, 0], [0, 0, 1]);
p = Polyhedra.polyhedron(v);
m = Polyhedra.Mesh(p);
X = mesh(m, color=:blue)
show(X)

using Oscar
RR, (x,y,z,l) = QQ["x","y","z"]
F0 = 37*x*y + 21*x*z + 22*y*z + 47*x + 2*y + 7*z + 35
F1 = 13*x^5*y^5*z^5 + 49*x + 44*y + 2*z + 14
p0 = newton_polytope(F0)
visualize(p0;VertexLabels = false)
p1 = newton_polytope(F1)
visualize(p0)

Polymake.visualize(p0.pm_polytope)