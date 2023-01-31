
needsPackage "tropicalToric"

verts0 = transpose matrix {{1,0,0},{0,1,0},{3,3,0},{0,0,1},{3,0,3},{0,3,3},{2,0,0},{0,2,0},{0,0,2}};
verts1 = transpose matrix {{3,0,0},{0,3,0},{5,5,5},{0,0,3},{2,0,0},{0,2,0},{0,0,2}};

verts0 = transpose matrix {{0,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,0},{0,1,1},{0,1,2}};
verts1 = transpose matrix {{1,0,0},{0,1,0},{2,2,2},{0,0,4},{0,2,0},{0,0,1},{0,0,0}};


verts0 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1}};
verts1 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,1},{0,0,1}};


needsPackage "tropicalToric"
P = minkowskiSum(convexHull verts0, convexHull verts1)
verts = matrix apply(entries vertices P, vert -> apply(vert, i->numerator i))
Y = makeSmooth normalToricVariety verts
mults = apply(0..length(rays Y)-1, i-> - min (entries ( matrix {(rays Y)_i} * verts ))_(0))

mults0 = apply(0..length(rays Y)-1, i-> - min (entries ( matrix {(rays Y)_i} * verts0 ))_(0))
c1E0 = sum toList apply(0..length(rays Y)-1, i->mults0_i*Y_i)
mults1 = apply(0..length(rays Y)-1, i-> - min (entries (matrix {(rays Y)_i} * verts1  ))_(0))
c1E1 = sum toList apply(0..length(rays Y)-1, i->mults1_i*Y_i)
(i, j, k) = toSequence positions(rays Y , r -> any({{1,0,0}, {0,1,0}, {0,0,1}}, e -> e == r)  );
c1H1 = Y_i ; c1H2 = Y_j ; c1H3 = Y_k ;
Wclass = c1H1 * c1H2 + c1H1 * c1H3 + c1H2 * c1H3 - (c1H1 + c1H2 + c1H3) * (c1E0 + c1E1) + c1E0 * c1E1 + c1E0 * c1E0 + c1E1 * c1E1 ;


allverts0 = fold( (v, w) -> v|w,  latticePoints convexHull verts0 )
allverts1 = fold( (v, w) -> v|w,  latticePoints convexHull verts1 )
R = QQ[X1,X2,X3];
f1 = sum apply(entries transpose allverts1, vert -> (random(500)*X1^(vert_0)*X2^(vert_1)*X3^(vert_2)));
f0 = sum apply(entries transpose allverts0, vert -> (random(500)*X1^(vert_0)*X2^(vert_1)*X3^(vert_2)));
IW = minors(2,jacobian(ideal(f0,f1)));
degCycle (c1E1 * Wclass)
degree (((ideal(f1, IW):ideal(X1)):ideal(X2)):ideal(X3))



A = ring Y;
newMon0 = vert -> random(200) * product toList apply(0..length(rayList)-1, i -> (A_i)^(mults0_i + (matrix({rayList_i})* vert)_0));
nverts0 = length(entries transpose verts0)
f0t = sum toList apply(0..nverts0 - 1, i -> newMon0 verts0_i )

newMon1 = vert -> random(200) * product toList apply(0..length(rayList)-1, i -> (A_i)^(mults1_i + (matrix({rayList_i})* vert)_0));
nverts1 = length(entries transpose verts1)
f1t = sum toList apply(0..nverts1 - 1, i -> newMon1 verts1_i )


inds = toList(set(1..length(rayList)-1) - set({i,j}))
Jt = jacobian(ideal(f0t,f1t));
J = submatrix'(Jt, inds,{})
IWt = minors(2,J);
codim ideal(f1t, IWt)




nVerts = rank target transpose  verts
rayList = rays(Y);
nRays = length(rayList);
vList = apply(0..nVerts -1 , i -> verts_i);
apply(0..nRays-1, i ->  sum toList select( vList  , vert ->  -1 * ( matrix {rayList_i} * vert)_0 == mults_i) )
 
mList = latticePoints convexHull verts0
O = apply((orbits Y)#1, orbit -> toSequence orbit)




F = diff(A_i, f0t) + diff(A_j, f0t) + diff(A_k, f0t) 
F1 = diff(A_1, f0t)
F0 = diff(A_0, f0t)
F4 = diff(A_4, f0t)

apply(O, orbit -> (orbit, apply((F0, F1, F4), pol -> sub(pol, {A_(orbit_0) => 0, A_(orbit_1) => 0}) == 0)))

someList =  select(O, orbit ->  sub(F4, {A_(orbit_0) => 0, A_(orbit_1) => 0}) == 0)
apply(someList, revealedFace)

select(O, orbit -> )


O_4
isAdmissible = (orbit, vert) -> all(orbit, rayIndex -> ( matrix{rayList_rayIndex} * vert_0)_0 == -mults0_rayIndex)
revealedFace = orbit -> select(mList,  vert -> isAdmissible(orbit, vert))
apply(O, orbit -> length revealedFace(orbit))

isAdmissiblePartial = (orbit,vert,ind) -> all(orbit,
(rayIndex -> if rayIndex == ind then ( matrix{rayList_rayIndex} * vert_0)_0 == 1
else ( matrix{rayList_rayIndex} * vert_0)_0 == -mults0_rayIndex))
revealedFacePartial = (ind, orbit) -> select(mList, vert -> isAdmissiblePartial(orbit, vert, ind))
apply(O, orbit -> length revealedFacePartial(i, orbit))









needsPackage "tropicalToric"
verts0 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1}};
verts1 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,1},{0,0,1}};
P = minkowskiSum(convexHull verts0, convexHull verts1)
verts = matrix apply(entries vertices P, vert -> apply(vert, i->numerator i))
Y = makeSmooth normalToricVariety verts


R = QQ[X1,X2,X3];
allverts0 = fold( (v, w) -> v|w,  latticePoints convexHull verts0 )
f0 = sum apply(entries transpose allverts0, vert -> (random(500)*X1^(vert_0)*X2^(vert_1)*X3^(vert_2)));

mults0 = apply(0..length(rays Y)-1, i-> - min (entries ( matrix {(rays Y)_i} * verts0 ))_(0))