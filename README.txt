Two implementation of closests point have been implemented.
The original kabsch code directly from the wikipedia which includes scale.
The second implementation is referred to as Julia.


This is a Julia implementation of a solution to a point-mapping problem:
given a reference set of points, and another set of points, how do you rigidly
rotate and translate the second set of points back onto the first, as nearly as possible?
Scaling the second set is not considered here.

This is sometimes called the Procrustes problem, or the Kabsch problem (when no translation is present).

My code is heavily indebted to Olga Sorkine's very readable explaination of using a
singular value decomposition (SVD) to solve the problem: [svd_rot.pdf] (http://igl.ethz.ch/projects/ARAP/svd_rot.pdf)

The Euler angle conventions are (I believe) consistent with much of the physics literature,
but are taken directly from
[Eric Weisstein's "Euler Angles" From MathWorld] (http://mathworld.wolfram.com/EulerAngles.html)

I also implemeted the algorithm in C++ using the Eigen library (http://eigen.tuxfamily.org/).

The C++ code can be compiled using a CMakeList.txt
The second implementation has been taken from align_scd.cpp code.
```g++ -I ../eigen/ align_svd.cpp -o align_svd```
but you may have to change the path to the Eigen libraries.

The C++ and Julia implementations give the same results for the test case I used (up to floating point error), but the actual affine transforms they use are different.
Not sure why this is the case.

There are series of tests in the tests directory. 

The transformation if from local geometry to world geometry.

The pitch,roll,yaw has been calculated for the transform. This not the pitch,roll and yaw of the system but the pitch,roll and between the systems. I am still checking to see if the results are correct
