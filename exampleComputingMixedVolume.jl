


using Oscar
R, (x,y,z,l) = QQ["x","y","z","l"]
f0 = 37*x*y + 21*x*z + 22*y*z + 47*x + 2*y + 7*z + 35
f1 = 13*x^5*y^5*z^5 + 49*x + 44*y + 2*z + 14
L = f0 + l * f1

polytopesLagrangian = newton_polytope.(jacobi_matrix(L))
polymakePols = [pol.pm_polytope for pol in polytopesLagrangian]
P1,P2,P3,P4 = polymakePols[1],polymakePols[2],polymakePols[3],polymakePols[4]
Polymake.polytope.mixed_volume(P1,P2,P3,P4)

function ∂(i,Δ::Polyhedron)   
    hyperplane = 
    intersect(hyperplane, Δ)
end
Cay = newton_polytope(L)
polytopes = [newton_polytope(f1), ∂(1,Cay), ∂(2,Cay), ∂(3,Cay)]
polymakePols = [polytopesLagrangian.pm_polytope for pol in polytopesLagrangian]
P1,P2,P3,P4 = polymakePols[1],polymakePols[2],polymakePols[3],polymakePols[4]
Polymake.polytope.mixed_volume(P1,P2,P3,P4)