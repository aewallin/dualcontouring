To test:
$ mkdir bld
$ cd bld
$ cmake ..
$ make 
$ ./dualcontour ../mechanic.dcf test.ply

Options:
--simplify 0.01  (octree simplification)
--nointer        (intersection-free algorithm)
--test           (run intersection tests after contouring)

This produces a test.ply file that can be viewed with meshlab.

original readme:
-------------------------------------
Dual Contouring Implementation in C++

Author: Tao Ju (with QEF code written by Scott Schaefer)
Updated: February 2011


I. What's included

/code		Source code and Microsoft Visualt Studio 6.0 project/workspace files
/data		A test file (mechanical part), in both .dcf and .ply formats


II. How to run

The dc.exe in the /code/release can be run by calling: 

>dc.exe mechanic.dcf out.ply

where out.ply stores the polygonal output.


III. File formats

The code can take in two kinds of input: .dcf (Dual Contouring Format) and 
.sog (Signed Octree with Geometry). Both formats store an octree grid with 
inside/outside signs. DCF contains intersection points and normals on 
grid edges, whereas SOG contains a single point location within each 
non-empty grid cell. Both formats can be produced from a polygonal 
model, via scan-conversion, using the Polymender software on my website:

http://www1.cse.wustl.edu/~taoju/code/polymender.htm

The detail formats are documented in the readme file of Polymender.


IV. Other notes.

Two algorithms are implemented in this code: the original dual contouring 
algorithm [Ju et al., Siggraph 2002] and the intersection-free 
extension [Ju et al., Pacific Graphics 2006]. You can switch between
 them in the main() function in dc.cpp. In addition, octree 
 simplification (guided by QEF errors) is also implemented, and can be 
 turned on in the main() function.

The use of all code is limited to non-profit research purposes only.
