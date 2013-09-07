/*

  Class for reading PLY (asc or binary) files
 * for PLY specification, see http://www.ics.uci.edu/~graphics/teaching/ply.html

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

#ifndef PLYREADER_H
#define PLYREADER_H

#include "GeoCommon.hpp"
#include "ModelReader.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

const char typeName[11][10] = {"char","uchar","short","ushort","int","uint","float","double", "uint8", "float32", "int32"} ;
const int typeSize[11] = {1,1,2,2,4,4,4,8,1,4,4} ;
const int totalTypes = 11 ;

class PLYReader : public ModelReader
{
public:
	// File handle
	FILE* fin ;

	// Number of triangles
	int numTrians, curtrian ;

	// Number of vertices
	int numVerts, curvert ;

	// Vertex array
	float* vertarray ;

	// Bounding box
	float min[3], max[3] ;
	float rawmin[3], rawmax[3] ;

	// Size of box
	float maxsize ;

	// Types and info
	int extraVertex ;

	// Temporary array for non-triangle faces
	int temparray[64] ;
	int tottri, curtri ;

	// Mode
	int mode ; // 0 for asc, 1 for little-endian, 2 for big endian

	// Mode
	int vmode ;

	// Offset for reading faces
	long offset ;

public:
	/// Constructor
	PLYReader( char* fname, float scl ) 
	{
		if ( ! ( fin = fopen( fname, "rb" ) ) )
		{
			printf("Unable to open file %s\n", fname) ;	
		};

		// Scan to get info
//		printf("Scanning file for bounding box...\n");
		
		// Parse header
		int vmode = 0, lines = 0 ;
		this->numTrians = this->numVerts = this->extraVertex = 0 ;
		char seps[] = " ,\t\n\r ";
		seps[5] = 10 ;
		while ( 1 )
		{
			char header[1024] ;
			fgets( header, 1024, fin ) ;
			// printf("%s\n", header) ;

			char* token = strtok( header, seps ) ;
			// token[ strlen(token) - 1 ] = '\0' ;
			// int comp = strcmp ( token, "end_header" ) ;
			// printf("%s %d\n", token, comp ) ;
			// char c = getchar() ;



			if ( ! strcmp ( token, "end_header" ) )
			{
				// Header finished
				// printf("End encountered!") ;
				break ;
			} else if ( ! strcmp( token, "format" ) ) 
			{
				// Format specification
				token = strtok ( NULL, seps ) ;
				if ( !strcmp ( token, "ascii" ) )
				{
					this->mode = 0 ;
				}
				else if ( ! strcmp ( token, "binary_big_endian") )
				{
					this->mode = 2 ;
				}
				else
				{
					this->mode = 1 ;
				}
				
			} else if ( ! strcmp( token, "element") )
			{
				vmode = 0 ;
				lines = 0 ;
				token = strtok( NULL, seps ) ;
				if ( !strcmp( token, "vertex" ) )
				{
					// Vertex count
					token = strtok( NULL, seps ) ;
					sscanf( token, "%d", &(this->numVerts) ) ;

					// Enter vertex mode
					vmode = 1 ;
				} else if ( !strcmp( token, "face" ) )
				{
					// Face count
					token = strtok( NULL, seps ) ;
					sscanf( token, "%d", &(this->numTrians) ) ;

					// Enter face mode
					vmode = 2 ;
				}
			} else if ( ! strcmp( token, "property" ) ) 
			{
				switch ( vmode )
				{
				case 1: // Vertex mode
					if ( lines >= 3 )
					{
						// Extra storage for each vertex
						token = strtok( NULL, seps ) ;

						for ( int i = 0 ; i < totalTypes ; i ++ )
						{
							if ( ! strcmp( token, typeName[i] ) )
							{
								this->extraVertex += typeSize[i] ;
								break ;
							}
						}
					}

					lines ++ ;
					break ;

				case 2: // Face mode

					break ;
				}
			}
		}
		//printf("Header info: %d vertices, %d faces, %d extra bytes per vertex. At %d\n", this->numVerts, this->numTrians, this->extraVertex, ftell( fin ) );
		this->vertarray = new float [ 3 * this->numVerts ] ;
		// char c = getchar() ;
		
		/*
		float *ok= new float[ numVerts * 3 ] ;
		char c = getchar() ;
		return ;
		
		this->vertarray = new float* [this->numVerts ] ;
		

		int i ;
		
		for ( i = 0 ; i < this->numVerts ; i ++ )
		{
			vertarray[i] = new float[4] ;
		}
		printf("Done");
		*/


		// Start reading vertices
		int i ;
		char temp[1024] ;
		if ( mode > 0 )
		{
			for ( i = 0 ; i < this->numVerts ; i ++ )
			{
				int haveread = 0 ;
				if ( ( haveread = fread( &(vertarray[ 3 * i ]), sizeof( float ), 3, fin ) ) != 3 )
				{

				  if( feof( fin ) )      
				  {
					 printf( "Read error %ld ", ftell( fin ) );
				  }

					printf("%d\n", haveread ) ;
					exit(0) ;
				}
				if ( mode == 2 )
				{
					// Flip bits
					for ( int j = 0 ; j < 3 ; j ++ )
					{
						flipBits32( & (vertarray[3*i + j]) ) ;
					}
				}
				/*
				printf("%f %f %f\n", vertarray[i][0], vertarray[i][1], vertarray[i][2] ) ;
				if ( i == 10 )
				{
					exit(0) ;
				}
				*/

				// Read extra data
				fread( temp, 1, this->extraVertex, fin ) ;
				if ( i == 0 )
				{
					for ( int j = 0 ; j < 3 ; j ++ )
					{
						rawmin[j] = rawmax[j] = vertarray[3 * i + j] ;
					}
				}
				else
				{
					for ( int k = 0 ; k < 3 ; k ++ )
					{
						if ( rawmin[k] > vertarray[3 * i + k] )
						{
							rawmin[k] = vertarray[3 * i + k] ;
						}
						if ( rawmax[k] < vertarray[3 * i + k] )
						{
							rawmax[k] = vertarray[3 * i + k] ;
						}
					}
				}

			}
		}
		else
		{
			// ASCII mode


			for ( i = 0 ; i < this->numVerts ; i ++ )
			{
				fgets ( temp, 1024, fin ) ;

				char seps[] = " ,\t\n";
				char* token = strtok( temp, seps ) ;

				sscanf( token, "%f", &(vertarray[ 3 * i ] ) ) ;
				token = strtok( NULL, seps ) ;
				sscanf( token, "%f", &(vertarray[ 3 * i + 1 ] ) ) ;
				token = strtok( NULL, seps ) ;
				sscanf( token, "%f", &(vertarray[ 3 * i + 2 ] ) ) ;

				if ( i == 0 )
				{
					for ( int j = 0 ; j < 3 ; j ++ )
					{
						rawmin[j] = rawmax[j] = vertarray[3 * i + j] ;
					}
				}
				else
				{
					for ( int k = 0 ; k < 3 ; k ++ )
					{
						if ( rawmin[k] > vertarray[3 * i + k] )
						{
							rawmin[k] = vertarray[3 * i + k] ;
						}
						if ( rawmax[k] < vertarray[3 * i + k] )
						{
							rawmax[k] = vertarray[3 * i + k] ;
						}
					}
				}

			}
		}


		maxsize = rawmax[0] - rawmin[0] ;
		if ( rawmax[1] - rawmin[1] > maxsize )
		{
			maxsize = rawmax[1] - rawmin[1] ;
		}
		if ( rawmax[2] - rawmin[2] > maxsize )
		{
			maxsize = rawmax[2] - rawmin[2] ;
		}

		for ( i = 0 ; i < 3 ; i ++ )
		{
			min[i] = ( rawmax[i] + rawmin[i] ) / 2 - maxsize / 2 ;
			max[i] = ( rawmax[i] + rawmin[i] ) / 2 + maxsize / 2 ;
		}

		// printf( "Low corner: %f %f %f   High corner %f %f %f \n", min[0], min[1], min[2], max[0], max[1], max[2] );
		// printf( "Maxdif: %f \n", maxsize );

		// Scale
		for ( i = 0 ; i < 3 ; i ++ )
		{
			min[i] -= maxsize * ( 1 / scl - 1 ) / 2 ;
		}
		maxsize *= ( 1 / scl ) ;

		// Reset 
		this->offset = ftell( fin ) ;
		curtrian = 0 ;

	};

	void reset( )
	{
		curtrian = 0 ;
		curvert = 0 ;

		curtri = 0 ;
		tottri = 0 ;

		fseek( fin, offset, SEEK_SET ) ;
	}

	/// Get next triangle
	Triangle* getNextTriangle( )
	{
		if ( curtrian == numTrians )
		{
			return NULL ;
		}	

		if ( curtri == tottri )
		{
			// Read next face
			unsigned char num ;


			if ( mode > 0 )
			{
				// Read number of vertices in the face
				fread( &num, sizeof( unsigned char ), 1, fin ) ;
				fread( temparray, sizeof( int ), num, fin ) ;
				if ( mode == 2 )
				{
					for ( int i = 0 ; i < num ; i ++ )
					{
						flipBits32( &(temparray[i]) ) ;
					}
				}
	//		printf("%d\n", tottri) ;
				// printf("\n") ;

				// if ( curtrian == 10 ) exit(1) ;
			}
			else
			{
				char temp[1024] ;
				fgets( temp, 1024, fin ) ;

				char seps[] = " ,\t\n";
				char* token = strtok( temp, seps ) ;

				sscanf( token, "%d", (int*)&num ) ;

				for ( int i = 0 ; i < num ; i ++ )
				{
					token = strtok( NULL, seps ) ;
					sscanf( token, "%d", &(temparray[i]) ) ;
				}
			}


			tottri = num - 2 ;
			curtri = 0 ;
			curtrian ++ ;
		}

		Triangle* t = new Triangle() ;
		for ( int i = 0 ; i < 3 ; i ++ )
		{
			t->vt[0][i] = vertarray[ 3 * temparray[0] + i ] ;
			t->vt[1][i] = vertarray[ 3 * temparray[curtri + 1] + i ] ;
			t->vt[2][i] = vertarray[ 3 * temparray[curtri + 2] + i ] ;
		}
		
		curtri ++ ;

		// printf("%d %d %d\n", temparray[0], temparray[curtri + 1], temparray[curtri + 2]);

		return t ;
	};
	
	/// Get next triangle
	int getNextTriangle( int ind[3] )
	{
		if ( curtrian == numTrians )
		{
			return 0 ;
		}		

		if ( curtri == tottri )
		{
			// Read next face
			unsigned char num ;

			if ( mode > 0 )
			{
				// Read number of vertices in the face
				fread( &num, sizeof( unsigned char ), 1, fin ) ;

				// Read array of indices
				fread( temparray, sizeof( int ), num, fin ) ;
				if ( mode == 2 )
				{
					for ( int i = 0 ; i < num ; i ++ )
					{
						flipBits32( &(temparray[i]) ) ;
						// printf("%d ", temparray[i]) ;
					}
				}
				// printf("\n") ;

				// if ( curtrian == 10 ) exit(1) ;
			}
			else
			{
				char temp[1024] ;
				fgets( temp, 1024, fin ) ;

				char seps[] = " ,\t\n";
				char* token = strtok( temp, seps ) ;

				sscanf( token, "%d", (int*)&num ) ;

				for ( int i = 0 ; i < num ; i ++ )
				{
					token = strtok( NULL, seps ) ;
					sscanf( token, "%d", &(temparray[i]) ) ;
				}
			}

			tottri = num - 2 ;
			curtri = 0 ;
			curtrian ++ ;
		}

		ind[0] = temparray[0] ;
		ind[1] = temparray[curtri + 1] ;
		ind[2] = temparray[curtri + 2] ;
		
		curtri ++ ;
		return 1 ;
	};

	void printInfo ( )
	{
		printf("Vertices: %d Polygons: %d\n",this->numVerts, this->numTrians) ;
	}


	/// Get bounding box
	float getBoundingBox ( float origin[3] )
	{
		for ( int i = 0 ; i < 3 ; i ++ )
		{
			origin[i] = min[i] ;
		}

		return maxsize ;
	};

	void getRawBoundingBox ( float low[3], float high[3] )
	{
		for ( int i = 0 ; i < 3 ; i ++ )
		{
			low[i] = rawmin[i] ;
			high[i] = rawmax[i] ;
		}
	}

	/// Get storage size
	int getMemory ( ) 
	{
		return sizeof ( class PLYReader ) + this->numVerts * sizeof( float ) * 3 ;	
	};

	// Flip bits

	void flipBits32 ( void *x )
	{
		unsigned char *temp = (unsigned char *)x;
		unsigned char swap;
		
		swap = temp [ 0 ];
		temp [ 0 ] = temp [ 3 ];
		temp [ 3 ] = swap;

		swap = temp [ 1 ];
		temp [ 1 ] = temp [ 2 ];
		temp [ 2 ] = swap;
	}

	void printBits32 ( void *x )
	{
		unsigned char *temp = (unsigned char *)x;

		printf("%2x %2x %2x %2x\n", temp[0], temp[1], temp[2], temp[3]) ;
	}
	
	int getNumTriangles()
	{
		return this->numTrians ;
	}

	int getNumVertices()
	{
		return this->numVerts ;
	}

	float* getVertex( int ind )
	{
		return &(this->vertarray[ 3 * ind ]) ;
	}

	void getNextVertex( float v[3] )
	{
		v[0] = vertarray[ 3 * curvert ]  ;
		v[1] = vertarray[ 3 * curvert + 1 ]  ;
		v[2] = vertarray[ 3 * curvert + 2 ]  ;

		curvert ++ ;
	}

	int getMode( )
	{
		return this->mode ;
	}

	void close()
	{
		fclose( this->fin ) ;
	}

};


#endif
