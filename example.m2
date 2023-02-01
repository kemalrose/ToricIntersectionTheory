

load "intersection.m2"

mons0 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1}};
mons1 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,1},{0,0,1}};

mons0 = transpose matrix {{0,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,0},{0,1,1},{0,1,2}};
mons1 = transpose matrix {{1,0,0},{0,1,0},{2,2,2},{0,0,4},{0,2,0},{0,0,1},{0,0,0}};

mons0 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{1,1,0},{0,0,1},{1,0,1},{0,1,1}};
mons1 = transpose matrix {{0,0,0},{1,0,0},{0,1,0},{5,5,5},{0,0,1}};

X = appropriateCompactification(mons0, mons1)
expectedNumberOfCriticalPoints(X,mons0,mons1)
numberOfCriticalPoints(mons0,mons1)

(c1E0,c1E1) = chernE(X,mons0,mons1);
(c1F1,c1F2, c1F3) = chernF(X,mons0,mons1);
WC = WClass(c1E0,c1E1,c1F1,c1F2, c1F3);
VC = c1E1;
degCycle (VC * WC)