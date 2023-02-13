


appropriateCompactification = (mons0, mons1) -> (
    P := minkowskiSum(convexHull mons0, convexHull mons1);
    newMons := matrix apply(entries vertices P, vert -> apply(vert, i->numerator i));
    Y := makeSmooth normalToricVariety newMons;
    Y
)

chernE = (Y,mons0,mons1) ->(
    mults0 := apply(0..length(rays Y)-1, i-> - min (entries ( matrix {(rays Y)_i} * mons0 ))_(0));
    c1E0 := sum toList apply(0..length(rays Y)-1, i->mults0_i*Y_i);
    mults1 := apply(0..length(rays Y)-1, i-> - min (entries (matrix {(rays Y)_i} * mons1  ))_(0));
    c1E1 := sum toList apply(0..length(rays Y)-1, i->mults1_i*Y_i);
    (c1E0,c1E1)
)

chernF = (Y,mons0,mons1) ->(
   (i, j, k) := toSequence positions(rays Y , r -> any({{1,0,0}, {0,1,0}, {0,0,1}}, e -> e == r)  );
    (Y_i, Y_j, Y_k)
)


WClass = (c1E0,c1E1,c1F1,c1F2, c1F3) ->(
    c1F1 * c1F2 + c1F1 * c1F3 + c1F2 * c1F3 - (c1F1 + c1F2 + c1F3) * (c1E0 + c1E1) + c1E0 * c1E1 + c1E0 * c1E0 + c1E1 * c1E1
)

expectedNumberOfCriticalPoints = (Y,mons0,mons1) ->(
    (c1E0,c1E1) := chernE(Y,mons0,mons1);
    (c1F1, c1F2, c1F3) := chernF(Y,mons0,mons1);
    WC := WClass(c1E0,c1E1,c1F1, c1F2, c1F3);
    degCycle (c1E1 * WC)
)


numberOfCriticalPoints = (mons0,mons1) -> (
    R = QQ[X1,X2,X3];
    f0 = sum apply(entries transpose mons0, mon -> (random(500)*X1^(mon_0)*X2^(mon_1)*X3^(mon_2)));
    f1 = sum apply(entries transpose mons1, mon -> (random(500)*X1^(mon_0)*X2^(mon_1)*X3^(mon_2)));
    IW = minors(2,jacobian(ideal(f0,f1)));
    (degree (((ideal(f1, IW):ideal(X1)):ideal(X2)):ideal(X3)), degree (ideal(f1, IW)))
)

fiberBundle = (Y,mons0,mons1) -> (
    
    normalToricVariety()
)
