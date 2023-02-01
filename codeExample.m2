
needsPackage "tropicalToric"
load "intersection.m2"

--This is the code accompanying example 5.14.
--mons0 and mons1 are the vertices of the poytopes P0 and P1 respectively.

mons0 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1}};
mons1 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{5,5,5},{0,0,1}};

-- X is an appropriate toric variety on which \EE and \FF are defined.

X = appropriateCompactification(mons0, mons1)

-- The following commands compute the chern classes of the line bundles L0, L1, H1, H2, H3 and the class of W

(c1E0,c1E1) = chernE(X,mons0,mons1);
(c1F1,c1F2, c1F3) = chernF(X,mons0,mons1);
WC = c1F1 * c1F2 + c1F1 * c1F3 + c1F2 * c1F3 - (c1F1 + c1F2 + c1F3) * (c1E0 + c1E1) + c1E0 * c1E1 + c1E0 * c1E0 + c1E1 * c1E1;

-- The following product of cycles is the number 69 
degCycle (c1E1 * WC)

-- This is equal to the number of critical points of a generic objective f0 and constraint f1 with fixed supports A0 and A1:

R = QQ[X1,X2,X3];
f0 = sum apply(entries transpose mons0, mon -> (random(500)*X1^(mon_0)*X2^(mon_1)*X3^(mon_2)));
f1 = sum apply(entries transpose mons1, mon -> (random(500)*X1^(mon_0)*X2^(mon_1)*X3^(mon_2)));
IW = minors(2,jacobian(ideal(f0,f1)));
(degree (((ideal(f1, IW):ideal(X1)):ideal(X2)):ideal(X3)), degree (ideal(f1, IW)))

