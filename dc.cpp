/*
  Copyright (C) 2011 Tao Ju

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "octree.hpp"
#include "PLYReader.hpp"
#include "PLYWriter.hpp"
#include "intersection.hpp"

#include <math.h>
#include <iostream>

//#define ALLOW_INTERSECTION

/*	Parameters
 *	argv[1]:	name of input file (.dcf format)
 *	argv[2]:	name of output file (.ply format)
 *	argv[3]:	(OPTIONAL) name of secondary output file (.ply format) 
 *              when using dual contouring, storing self-intersecting triangles.
*/

int main( int args, char* argv[] )
{
	// Reading input file
	std::cout << " input file: " << argv[1] << "\n";
	Octree* mytree = new Octree( argv[1] ) ;

	// Octree simplification; feel free to change the simplification threshold.
	// mytree->simplify( .001f ) ;

#ifdef ALLOW_INTERSECTION

	// Dual contouring [Ju et al. 2002]
	mytree->genContour( argv[2] ) ;

	// Pairwise intersection test - may take a while
	// int num = Intersection::testIntersection( argv[2], argv[3] ) ;
	// printf("%d intersections found!\n", num) ;

#else
	std::cout << " Intersection-free algorithm! \n";

	// Intersection-free dual contouring [Ju et al. 2006]
	mytree->genContourNoInter2( argv[2] ) ;

#endif
}

