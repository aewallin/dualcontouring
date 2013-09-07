/*

  Virtual class for input file readers

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



#ifndef MODELREADER_H
#define MODELREADER_H

#include "GeoCommon.hpp"

class ModelReader
{
public:
	/// Constructor
	ModelReader(){} ;
	/// Get next triangle
	virtual Triangle* getNextTriangle( ) = 0 ;
	virtual int getNextTriangle( int t[3] ) = 0 ;
	/// Get bounding box
	virtual float getBoundingBox ( float origin[3] ) = 0 ;
	/// Get number of triangles
	virtual int getNumTriangles ( ) = 0 ;
	/// Get storage size
	virtual int getMemory ( ) = 0 ;
	/// Reset file reading location
	virtual void reset( ) = 0 ;
	/// For explicit vertex models
	virtual int getNumVertices( ) = 0 ;
	virtual void getNextVertex( float v[3] ) = 0 ;
	virtual void printInfo ( ) = 0 ;
};


#endif
