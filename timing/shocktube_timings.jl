# Timings for BATSRUS explicit shocktube tests.
#
# Platform:
# Yuxi Linux Redhat
# Stampede2
# Frontera
#
# Hongyang Zhou, hyzhou@umich.edu 07/01/2019

using PyPlot

# noopenmp: 1,2,4 mpi
# openmp:   (1,1) (1,2) (2,1) (1,4) (2,2) (4,1)
# Yuxi Linux ifort 2018, Yuxi Linux gfortran 4.3,
# Stampede2 ifort, Stampede2 gfortran,
# Frontera ifort, Frontera gfortran
times = [0.69 0.37 0.21 1.00 0.56 0.56 0.31 0.30 0.29;
        0.82 0.43 0.24 0.85 0.46 0.44 0.26 0.25 0.25;
        0.79 0.48 0.30 1.14 1.10 0.66 0.68 0.49 0.44;
        1.87 0.55 0.38 1.88 1.11 0.57 0.69 0.70 0.37;
        0.63 0.36 0.20 0.94 0.91 0.51 0.57 0.36 0.28;
        0.62 0.34 0.21 0.63 0.36 0.38 0.22 0.21 0.21]

proc = [1; 2]

plt.figure(figsize=(5,5))
plt.plot(proc,times[:,[1,4]]',"o",Markersize=5,alpha=0.8)

plt.legend(loc="best",
   ("Redhat ifort19","Redhat gfortran4.8",
   "Stampede2 ifort18","Stampede2 gfortran7.3",
   "Frontera ifort19","Frontera gfortran9.0") )
plt.grid(true)
plt.title("BATSRUS Shocktube Test Performance, nProc*nthread=1")
plt.xlabel("mpi process")
plt.xticks(proc,("without OpenMP","with OpenMP"))
plt.ylabel("timing [s]")


proc = [1;2;3]
plt.figure(figsize=(5,5))
plt.plot(proc,times[:,[2,5,6]]',"-o",Markersize=5,alpha=0.8)

plt.legend(loc="best",
   ("Redhat ifort19","Redhat gfortran4.8",
   "Stampede2 ifort18","Stampede2 gfortran7.3",
   "Frontera ifort19","Frontera gfortran9.0"), framealpha=0.3)
plt.grid(true)
plt.title("BATSRUS Shocktube Test Performance, nProc*nthread=2")
plt.xticks(proc,("without OpenMP","with OpenMP (1,2)","with OpenMP (2,1)"))
plt.ylabel("timing [s]")


proc = [1;2;3;4]
plt.figure(figsize=(8,5))
plt.plot(proc,times[:,[3,9,8,7]]',"-o",Markersize=5,alpha=0.8)

plt.legend(loc="best",
   ("Redhat ifort19","Redhat gfortran4.8",
   "Stampede2 ifort18","Stampede2 gfortran7.3",
   "Frontera ifort19","Frontera gfortran9.0"), framealpha=0.5 )
plt.grid(true)
plt.title("BATSRUS Shocktube Test Performance, nProc*nthread=4")
plt.xticks(proc,("without OpenMP","with OpenMP (4,1)","with OpenMP (2,2)",
   "with OpenMP (1,4)"))
plt.ylabel("timing [s]")
