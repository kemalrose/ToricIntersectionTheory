

using Oscar

P0 = convex_hull([[0,0,0],[2,0,0],[0,2,0]])
P1 = convex_hull([[2,0,2],[0,2,2],[2,2,2]])
#Note: P1 - D_{e_1} is virtual.
Cay = convex_hull(P0,P1)
visualize(Cay)


facet_list = facets(Cay)
A, b = halfspace_matrix_pair(facet_list)
#construct new polyhedron by intersecting with halspace {x : xi ≥ 1}
halfspace = Polyhedron([-1,0,0],-1)
facets(halfspace)
∂1Cay = intersect(Cay, halfspace)
visualize(∂1Cay)

#Now compute the rays and their multiplicities for both Cay and ∂1Cay.
A1, b1 = halfspace_matrix_pair(facets(∂1Cay))
A2, b2 = halfspace_matrix_pair(facets(Cay))

#Note: up to permutation A1 and A2 coincide. b1 and b2 differ only in the single entry corresponding to the vector e1.


