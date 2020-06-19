# voronoi
Ising spin models on voronoi tesselation

build in the usual cmake way after installing morphologica

run using e.g., ./build/voronoi.cpp 10 5

which will run the ising model on a tesselation created from 10 random seed points, at 5 different temperature values (in the range 0 to 50 by default).

Then run ising.py to see results (energy etc. versus temperature) and voronoi.py which will show the tesselation (for the final sim run), with an overlay of info from one cell (in cell.h5) which I was using for debugging.
