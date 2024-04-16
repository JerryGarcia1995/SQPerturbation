1. Compile LSCoupling.cpp with "Eigen" library: e.g. mpicxx LSCoupling.cpp -o LSCoupling.x -std=c++11 -I/Eigen
2. Run LSCoupling.x. A good example is given in oper_LSCoupling.x. Customize this to your architecture.
3. oper_LSCoupling.x will give many cases of eigenvalues/eigenvectors for f^1, ..., f^13.
   The recommendations are given below (to fit best to the Dieke diagram)
   f01Ce3+_JH0.1.dat,
   f02Pr3+_JH0.6.dat,
   f03Nd3+_JH0.7.dat,
   f04Pm3+_JH0.8.dat,
   f05Sm3+_JH0.7.dat,
   f06Eu3+_JH0.9.dat,
   f07Gd3+_JH1.0.dat,
   f08Tb3+_JH1.0.dat,
   f09Dy3+_JH1.1.dat,
   f10Ho3+_JH0.8.dat,
   f11Er3+_JH0.9.dat,
   f12Tm3+_JH0.5.dat, and
   f13Yb3+_JH0.1.dat
4. Compile SQPerturbation.cpp: e.g. mpicxx SQPerturbation.cpp -o SQPerturbation.x -std=c++11
5. Run SQPerturbation.x. Good examples are given in oper_f****3+G*.sh.
   For example, oper_f01Ce3+G7.sh will perform second perturbation calculations for Gamma_7 of f^1: f^1 f^1 -> f^0 f^2 -> f^1 f^1
   Note that HtSH_Oh.pfpi.-0.7.txt provides the Slater-Koster parameters for pfpi/pfsigma = -0.7.
   Note that ManualMPI_f *.tsv is supposed to maximize the parallel efficiency for Ohtaka (supercomputer at ISSP). If you are not sure of this, set ManualMPI to 1 and ManualMPIFilename to any. (If ManualMPI = 1, ManualMPIFilename will be just passed by.)
7. The output files for SQPerturbation.x will show the coupling constants. Meanwhile, you can also analyze which intermediate states contribute to the coupling constants most.
   Compile SQPerturbationDecompositionAnalyzer.cpp: e.g. mpicxx SQPerturbationDecompositionAnalyzer.cpp -o SQPerturbationDecompositionAnalyzer.x -std=c++11

Note that the final results are given in results_raw_data.zip.

I apologize for the short explanation. But if you are a geek, you will have no problem, I believe.
