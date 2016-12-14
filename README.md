MandelbrotFortran
=================

###Introduction
---

This project creates an image of a [Mandelbrot fractal](http://en.wikipedia.org/wiki/Mandelbrot_set). This repository contains 
three separate codes: a serial version, an OpenMP version, and a hybrid OpenMP and MPI version.

This project uses the [Distance Estimator](http://mrob.com/pub/muency/distanceestimator.html) method to calculate the Mandelbrot 
set. Each pixel requires its own independent unit of work. Because these units of work are independent, it is a prime candidate 
for parallelizaion. 

The code will perform an iterative calculation at each point that will determine weather if the point 
[escapes the set](http://en.wikipedia.org/wiki/Mandelbrot_set#Formal_definition). If the point doesn't escape, then it is 
considered part of the Mandelbrot set and the pixel is painted black.

The Distance Estimator comes into play if the point escapes very slowly (takes many iterations before it escapes). For 
visualization purposes, we consider these points to be part of the set. By doing this, it will "reveal" more of the set than the 
standard method. 

####Serial
This code is located in `serial_mandelbrot.f90`. At the beginning of the serial code, a 2 dimensional grid is allocated in memory. 
This grid will hold either a 1 or 0 in each element to represent whether a point is in the set or not. 

The code then iterates over each point in the grid and performs the Distance Estimator calculation. When the calculation determines 
if the point escapes, it breaks out of the loop and writes a "0" in the grid for that point. Else, it will infer that the point does 
not escape and will write a "1" in the grid. 

After iterating over each point, the grid is then passed to a function that will write out a pgm image to the current working 
directory. The function is code by John Burkardt and can be obtained [here](http://people.sc.fsu.edu/~jburkardt/f_src/pgma_io/pgma_io.html)

####Compiling
To compile the serial code use your fortran compiler, for example

  gfortran serial_mandelbrot.f90 -o serial_mandelbrot

The code can then be run

  ./serial_mandelbrot

and once it completes should put a picture file in the running directory. The picture can be opened in many viewers, such as 
[GIMP](https://www.gimp.org/) or [Gwenview](https://userbase.kde.org/Gwenview)
