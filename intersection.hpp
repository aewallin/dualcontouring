/*

  Copyright (C) 2011  Tao Ju

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

#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "GeoCommon.hpp"
#include "PLYReader.hpp"
#include "PLYWriter.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Intersection
{
public:
	Intersection(){} ;
     // c = a x b
	static void cross( float a[3], float b[3], float c[3] )
	{
		c[0] = a[1] * b[2] - a[2] * b[1] ;
		c[1] = a[2] * b[0] - a[0] * b[2] ;
		c[2] = a[0] * b[1] - a[1] * b[0] ;
	}
	// a dot b
	static float dot( float a[3], float b[3] )
	{
		return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] ;
	}
	
	
	static int separating( float axes[3], Triangle* t1, Triangle* t2 )
	{
		//float mag1 = sqrt(dot( axes, axes )) ;
		//axes[0] /= mag1 ;
		//axes[1] /= mag1 ;
		//axes[2] /= mag1 ;
		
		float min1 = dot( axes, t1->vt[0] ), max1;
		max1 = min1 ;
		float min2 = dot( axes, t2->vt[0] ), max2;
		max2 = min2 ;
		
		float temp ;
		for ( int i = 1 ; i < 3 ; i ++ )
		{
			temp = dot( axes, t1->vt[i] ) ;
			if ( temp < min1 )
			{
				min1 = temp ;
			}
			else if ( temp > max1 )
			{
				max1 = temp ;
			}
			
			temp = dot( axes, t2->vt[i] ) ;
			if ( temp < min2 )
			{
				min2 = temp ;
			}
			else if ( temp > max2 )
			{
				max2 = temp ;
			}
		}
		
		if ( min1 >= max2 || min2 >= max1 )
		{
			return 1 ;
		}
		else
		{
			return 0 ;
		}
	}

	static int separating( float axes[3], Triangle* t1, Triangle* t2, int printout )
	{
		//float mag1 = sqrt(dot( axes, axes )) ;
		//axes[0] /= mag1 ;
		//axes[1] /= mag1 ;
		//axes[2] /= mag1 ;
		
		float min1 = dot( axes, t1->vt[0] ), max1;
		max1 = min1 ;
		float min2 = dot( axes, t2->vt[0] ), max2;
		max2 = min2 ;
		
		float temp ;
		for ( int i = 1 ; i < 3 ; i ++ )
		{
			temp = dot( axes, t1->vt[i] ) ;
			if ( temp < min1 )
			{
				min1 = temp ;
			}
			else if ( temp > max1 )
			{
				max1 = temp ;
			}
			
			temp = dot( axes, t2->vt[i] ) ;
			if ( temp < min2 )
			{
				min2 = temp ;
			}
			else if ( temp > max2 )
			{
				max2 = temp ;
			}
		}

		if ( printout )
		{
					for ( int k = 0 ; k < 3 ; k ++ )
					{
						printf("(%f %f %f),", t1->vt[k][0], t1->vt[k][1],t1->vt[k][2]);
					}
					printf("\n");
					for ( int k = 0 ; k < 3 ; k ++ )
					{
						printf("(%f %f %f),", t2->vt[k][0], t2->vt[k][1],t2->vt[k][2]);
					}

			printf("\n(%f, %f), (%f, %f)\n", min1, max1, min2, max2) ;
		}
		
		if ( min1 >= max2 || min2 >= max1 )
		{
			return 1 ;
		}
		else
		{
			return 0 ;
		}
	}
	
	static int separating( float axes[3], Triangle* t1, float v1[3], float v2[3] )
	{
		float mag1 = sqrt(dot( axes, axes )) ;
		axes[0] /= mag1 ;
		axes[1] /= mag1 ;
		axes[2] /= mag1 ;
		
		float min1 = dot( axes, t1->vt[0] ), max1 = min1 ;
		float temp ;
		for ( int i = 1 ; i < 3 ; i ++ )
		{
			temp = dot( axes, t1->vt[i] ) ;
			if ( temp < min1 )
			{
				min1 = temp ;
			}
			else if ( temp > max1 )
			{
				max1 = temp ;
			}
		}

		float min2 = dot( axes, v1 ) ;
		float max2 = dot( axes, v2 ) ;
		if ( min2 > max2 )
		{
			temp = min2 ;
			min2 = max2 ; 
			max2 = temp ;
		}

		
		if ( min1 > max2 + 0.00001f || min2 > max1 + 0.00001f )
		{
			return 1 ;
		}
		else
		{
			return 0 ;
		}
	}
	
	static int testIntersection( Triangle* t1, Triangle* t2 )
	{
		// Two face normals
		float v1[3][3] = {{ t1->vt[1][0] - t1->vt[0][0], t1->vt[1][1] - t1->vt[0][1], t1->vt[1][2] - t1->vt[0][2] },
		{ t1->vt[2][0] - t1->vt[1][0], t1->vt[2][1] - t1->vt[1][1], t1->vt[2][2] - t1->vt[1][2] },
		{ t1->vt[0][0] - t1->vt[2][0], t1->vt[0][1] - t1->vt[2][1], t1->vt[0][2] - t1->vt[2][2] }};
		float v2[3][3] = {{ t2->vt[1][0] - t2->vt[0][0], t2->vt[1][1] - t2->vt[0][1], t2->vt[1][2] - t2->vt[0][2] },
		{ t2->vt[2][0] - t2->vt[1][0], t2->vt[2][1] - t2->vt[1][1], t2->vt[2][2] - t2->vt[1][2] },
		{ t2->vt[0][0] - t2->vt[2][0], t2->vt[0][1] - t2->vt[2][1], t2->vt[0][2] - t2->vt[2][2] } };
		float n1[3], n2[3] ; 
		cross( v1[0], v1[1], n1 ) ;
		cross( v2[0], v2[1], n2 ) ;
		
		float n[3] ;
		cross( n1, n2, n ) ;
		int i, j ;
		if ( n[0] == 0 && n[1] == 0 && n[2] == 0 )
		{
			// Co-planar
			
			if ( dot( n1, t1->vt[0] ) != dot( n1, t2->vt[0] ) )
			{
				return 0 ;
			}
			
			float axes[3] ;
			for ( i = 0 ; i < 3 ; i ++ )
			{
				cross( n1, v1[i], axes ) ;
				if ( separating( axes, t1, t2 ) )
				{
					return 0 ;
				}
			}
			for ( i = 0 ; i < 3 ; i ++ )
			{
				cross( n2, v2[i], axes ) ;
				if ( separating( axes, t1, t2 ) )
				{
					return 0 ;
				}
			}
		}
		else
		{
			// Non co-planar
			if ( separating( n1, t1, t2 ) || separating( n2, t1, t2 ) )
			{
				return 0 ;
			}
			
			float axes[3] ;
			for ( i = 0 ; i < 3 ; i ++ )
				for ( j = 0 ; j < 3 ; j ++ )
				{
					cross( v1[i], v2[j], axes ) ;
					if ( separating( axes, t1, t2 ) )
					{
						return 0 ;
					}
				}
		}
		
		return 1 ;
	}

	static int testIntersection( Triangle* t1, BoundingBox b1, Triangle* t2, BoundingBox b2, int printout )
	{
		// Bounding box
		if ( b1.begin.x > b2.end.x || b1.begin.y > b2.end.y || b1.begin.z > b2.end.z ||
			b2.begin.x > b1.end.x || b2.begin.y > b1.end.y || b2.begin.z > b1.end.z )
		{
			return 0 ;
		}

		// Two face normals
		float v1[3][3] = {{ t1->vt[1][0] - t1->vt[0][0], t1->vt[1][1] - t1->vt[0][1], t1->vt[1][2] - t1->vt[0][2] },
		{ t1->vt[2][0] - t1->vt[1][0], t1->vt[2][1] - t1->vt[1][1], t1->vt[2][2] - t1->vt[1][2] },
		{ t1->vt[0][0] - t1->vt[2][0], t1->vt[0][1] - t1->vt[2][1], t1->vt[0][2] - t1->vt[2][2] }};
		float v2[3][3] = {{ t2->vt[1][0] - t2->vt[0][0], t2->vt[1][1] - t2->vt[0][1], t2->vt[1][2] - t2->vt[0][2] },
		{ t2->vt[2][0] - t2->vt[1][0], t2->vt[2][1] - t2->vt[1][1], t2->vt[2][2] - t2->vt[1][2] },
		{ t2->vt[0][0] - t2->vt[2][0], t2->vt[0][1] - t2->vt[2][1], t2->vt[0][2] - t2->vt[2][2] } };
		float n1[3], n2[3] ; 
		cross( v1[0], v1[1], n1 ) ;
		cross( v2[0], v2[1], n2 ) ;
		
		float n[3] ;
		cross( n1, n2, n ) ;
		int i, j ;
		if ( n[0] == 0 && n[1] == 0 && n[2] == 0 )
		{
			// Co-planar, regard it as not intersecting -- Hack!
			return 0 ;

			if ( printout )
			{
				printf("Co-planar!\n") ; 
			}
			if ( dot( n1, t1->vt[0] ) != dot( n1, t2->vt[0] ) )
			{
				return 0 ;
			}
			
			float axes[3] ;
			for ( i = 0 ; i < 3 ; i ++ )
			{
				cross( n1, v1[i], axes ) ;
				if ( separating( axes, t1, t2, printout ) )
				{
					return 0 ;
				}
			}
			for ( i = 0 ; i < 3 ; i ++ )
			{
				cross( n2, v2[i], axes ) ;
				if ( separating( axes, t1, t2, printout ) )
				{
					return 0 ;
				}
			}
		}
		else
		{
			// Non co-planar
			if ( separating( n1, t1, t2, printout ) || separating( n2, t1, t2, printout ) )
			{
				return 0 ;
			}
			
			float axes[3] ;
			for ( i = 0 ; i < 3 ; i ++ )
				for ( j = 0 ; j < 3 ; j ++ )
				{
					cross( v1[i], v2[j], axes ) ;
					if ( separating( axes, t1, t2, printout ) )
					{
						return 0 ;
					}
				}
		}
		
		return 1 ;
	}
	
	static int testIntersection( char* fname, char* outname )
	{
		Triangle *intertris[100000];
		
		// Read triangles
		PLYReader* myreader = new PLYReader( fname, 1 ) ;
		int numpoly = myreader->getNumTriangles() ;
		int num = 0 ;
		Triangle** tris = new Triangle *[ numpoly ] ;
		Triangle* tri = NULL ;
		myreader->reset() ;
		while( (tri = myreader->getNextTriangle()) != NULL )
		{
			tris[ num ] = tri ;
			num ++ ;
		}
		myreader->close() ;
		printf("Reading %d polygons, %d triangles.\n", numpoly, num ) ;
		
		// Build bounding boxes
		printf("Building bounding boxes...\n") ;
		int i ;
		BoundingBox* boxes = new BoundingBox[ num ] ;
		for ( i = 0 ; i < num ; i ++ )
		{
			Triangle* t = tris[ i ] ;
			
			boxes[i].begin.x = 100000 ;
			boxes[i].begin.y = 100000 ;
			boxes[i].begin.z = 100000 ;
			boxes[i].end.x = -100000 ;
			boxes[i].end.y = -100000 ;
			boxes[i].end.z = -100000 ;
			
			for ( int j = 0 ; j < 3 ; j ++ )
			{
				if ( t->vt[j][0] < boxes[i].begin.x )
				{
					boxes[i].begin.x = t->vt[j][0] ;
				}
				if ( t->vt[j][0] > boxes[i].end.x )
				{
					boxes[i].end.x = t->vt[j][0] ;
				}
				
				if ( t->vt[j][1] < boxes[i].begin.y )
				{
					boxes[i].begin.y = t->vt[j][1] ;
				}
				if ( t->vt[j][1] > boxes[i].end.y )
				{
					boxes[i].end.y = t->vt[j][1] ;
				}
				
				if ( t->vt[j][2] < boxes[i].begin.z )
				{
					boxes[i].begin.z = t->vt[j][2] ;
				}
				if ( t->vt[j][2] > boxes[i].end.z )
				{
					boxes[i].end.z = t->vt[j][2] ;
				}
			}
			
		}
		
		// Test intersections
		printf("Pairwise searching...\n") ;
		int numinters = 0 ;
		for ( i = 0 ; i < num ; i ++ )
		{
			for ( int j = 0 ; j < i ; j ++ )
			{
				if ( testIntersection( tris[i], boxes[i], tris[j], boxes[j], 0 ) )
				{
					intertris[ 2 * numinters ] = tris[i] ;
					intertris[ 2 * numinters + 1 ] = tris[j] ;
					numinters ++ ;
				
					/*
					printf("First tri: ");
					for ( int k = 0 ; k < 3 ; k ++ )
					{
						printf("(%f %f %f),", tris[i]->vt[k][0], tris[i]->vt[k][1],tris[i]->vt[k][2]);
					}
					printf("\nSecond tri: ");
					for ( k = 0 ; k < 3 ; k ++ )
					{
						printf("(%f %f %f),", tris[j]->vt[k][0], tris[j]->vt[k][1],tris[j]->vt[k][2]);
					}
					// printf("\n") ;
					*/
				}
			}
		}
		
		// Write out
		printf("Write out...\n") ;
		if ( numinters > 0  )
		{
			FILE* intout = fopen( outname, "wb" ) ;
			PLYWriter::writeHeader( intout, 6 * numinters, 2 * numinters ) ;
			Triangle* t ;
			
			// Write vertices
			for ( int i = 0 ; i < numinters ; i ++ )
			{
				t = intertris[ 2 * i ] ;
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					PLYWriter::writeVertex( intout, t->vt[j] ) ;
				}
				t = intertris[ 2 * i + 1 ] ;
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					PLYWriter::writeVertex( intout, t->vt[j] ) ;
				}
			}
			
			// Write triangles
			for ( int i = 0 ; i < numinters ; i ++ )
			{
				int tind1[] = { 6 * i, 6 * i + 1, 6 * i + 2 } ;
				PLYWriter::writeFace( intout, 3, tind1 ) ;
				int tind2[] = { 6 * i + 3, 6 * i + 4, 6 * i + 5 } ;
				PLYWriter::writeFace( intout, 3, tind2 ) ;
			}
			
			fclose( intout ) ;
		}
		
		return numinters ;
	}

};


#endif
