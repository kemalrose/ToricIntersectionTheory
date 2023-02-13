


using Oscar

A,(z1,z2,l) = PolynomialRing(QQ,["z1","z2","l"])
function random_pol(verts)
    allverts = Vector{Int64}.(lattice_points(convex_hull(verts)))
    f = rand(-500:500, length(allverts))'*[z1^pt[1]*z2^pt[2] for pt in allverts]
    f
end

function nr_crits(verts0,verts1)
    f0, f1 = random_pol.([verts0, verts1])
    L =  f0 + l * f1
    ∇L = jacobi_matrix(L)
    cIdeal = ideal(A,∇L.entries[:])
    vdim(quo(A,cIdeal)[1])
end

function scaled_nr_crits(verts0, verts1, c)
    f0, f1 = random_pol.([c*verts0, c*verts1])
    L = f0 + l * f1
    ∇L = [L,L,L]
    for i in 1:c
        ∇L = [derivative(∇L[j],j) for j in 1:length(∇L)]
    end
    cIdeal = ideal(A,∇L)
    vdim(quo(A,cIdeal)[1])
end

function ∂(i, pol)
    A = [0,0,0]
    A[i] = -1
    b = -1
    half_plane =Polyhedron(A,b)
    intersect(halp_plane, pol)
end



verts0 = [[0,0],[1,0],[0,1],[1,4]]
verts1 = [[0,0],[2,0],[0,2],[4,3]]

f0, f1 = random_pol.([verts0, verts1])
L =  f0 + l * f1
∇L = jacobi_matrix(L)
cIdeal = ideal(A,∇L.entries[:])
vdim(quo(A,cIdeal)[1])
polytopes = newton_polytope.(jacobi_matrix(L))
polymakePols = [polytopes[i].pm_polytope for i = 1:3]
P1,P2,P3 = polymakePols[1],polymakePols[2],polymakePols[3]
Polymake.polytope.mixed_volume(polymakePols...)

Cay = newton_polytope(L)
∇L = jacobi_matrix(L)
[∂(i, pol)]



V = ideal(A,[f1])
J = jacobi_matrix([f0,f1])[1:2,:]
W = ideal(A, minors(J,2) )
VcapW = V+W 
vdim(quo(A,V+W)[1])



using Oscar
R, (x,y,z,l) = QQ["x","y","z","l"]
f0 = 37*x*y + 21*x*z + 22*y*z + 47*x + 2*y + 7*z + 35
f1 = 13*x^5*y^5*z^5 + 49*x + 44*y + 2*z + 14
L = f0 + l * f1

polytopes = newton_polytope.(jacobi_matrix(L))
polymakePols = [polytopes[i].pm_polytope for i = 1:4]
P1,P2,P3,P4 = polymakePols[1],polymakePols[2],polymakePols[3],polymakePols[4]
Polymake.polytope.mixed_volume(polymakePols...)
