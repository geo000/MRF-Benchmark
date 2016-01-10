READ ME

To use the Dual Decomposition algorithm, just use the executable as usual, but with the option "-a 6".
Note : "-a 7" will execute all algorithms, including the dual decomposition.
The code should be compiled if "make" is typed in the MRF-benchmark folder. However, since we had to reinstall some part of imageLib, it may not happen. The executables will work, though.

E.G. : 
./mrfstereo -a 6 data/tsukuba-imL.png  data/tsukuba-imR.png tsukuba-DD.png
