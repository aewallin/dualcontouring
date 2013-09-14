/*

  Implementations of Octree member functions.

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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "octree.hpp"
#include "PLYWriter.hpp"

Octree::Octree( char* fname )
{
	// Recognize file format
	if ( strstr( fname, ".sog" ) != NULL || strstr( fname, ".SOG" ) != NULL ) {
		printf("Reading SOG file format.\n") ;
		this->hasQEF = 0 ;
		readSOG( fname ) ;
	} 
	else if ( strstr( fname, ".dcf" ) != NULL || strstr( fname, ".DCF" ) != NULL ) {
		printf("Reading DCF file format.\n") ;
		this->hasQEF = 1 ;
		readDCF( fname ) ;
	}
	else {
		printf("Wrong input format. Must be SOG/DCF.\n");
		exit(0) ;
	}

}

void Octree::simplify( float thresh ) {
	if ( this->hasQEF ) {
		int st[3] = {0,0,0} ;
		this->root = simplify( this->root, st, this->dimen, thresh ) ;
	}
}

// simplify by collapsing nodes where the parent node QEF solution is good enough
OctreeNode* Octree::simplify( OctreeNode* node, int st[3], int len, float thresh ) {
	if ( node == NULL )
		return NULL ;

	NodeType type = node->getType() ;

	if ( type == INTERNAL ) { 	// Internal node
		InternalNode* inode = (InternalNode*)node ;
		int simple = 1 ;

		float ata[6] = { 0, 0, 0, 0, 0, 0 };
		float atb[3] = { 0, 0, 0 } ;
		float pt[3] = { 0, 0, 0 } ;
		float mp[3] = { 0, 0, 0 } ;
		float btb = 0 ;
		int signs[8] = {-1,-1,-1,-1,-1,-1,-1,-1} ;
		int midsign = -1 ;

		int nlen = len / 2 ;
		int nst[3] ;
		int ec = 0 ;
		int ht ;

		for ( int i = 0 ; i < 8 ; i ++ ) {
			nst[0] = st[0] + vertMap[i][0] * nlen ;
			nst[1] = st[1] + vertMap[i][1] * nlen ;
			nst[2] = st[2] + vertMap[i][2] * nlen ;

			inode->child[i] = simplify( inode->child[i], nst, nlen, thresh ) ;
			
			if ( inode->child[i] != NULL ) {
				if ( inode->child[i]->getType() == 0 ) {
					simple = 0 ;
				}
				else if ( inode->child[i]->getType() == 1 ) {
					LeafNode* lnode = (LeafNode *) inode->child[i] ;
					ht = lnode->height ;

					for ( int j = 0 ; j < 6 ; j ++ )
						ata[j] += lnode->ata[j] ;

					for ( int j = 0 ; j < 3 ; j ++ ) {
						atb[j] += lnode->atb[j] ;
						pt[j] += lnode->mp[j] ;
					}
					if ( lnode->mp[0] == 0 )
						printf("%f %f %f, Height: %d\n", lnode->mp[0], lnode->mp[1], lnode->mp[2], ht) ;

					btb += lnode->btb ;
					ec ++ ;

					midsign = lnode->getSign( 7 - i ) ;
					signs[i] = lnode->getSign( i ) ;
				}
				else {
					PseudoLeafNode* pnode = (PseudoLeafNode *) inode->child[i] ;
					ht = pnode->height ;

					for ( int j = 0 ; j < 6 ; j ++ )
						ata[j] += pnode->ata[j] ;

					for ( int j = 0 ; j < 3 ; j ++ ) {
						atb[j] += pnode->atb[j] ;
						pt[j] += pnode->mp[j] ;
					}
					btb += pnode->btb ;
					ec ++ ;

					midsign = pnode->getSign( 7 - i ) ;
					signs[i] = pnode->getSign( i ) ;
				}
			}
		}

		if ( simple )
		{
			// Ok, let's collapse
			if ( ec == 0 ) {
				delete node ;
				return NULL ;
			}
			else {
				pt[0] = pt[0] / ec ;
				pt[1] = pt[1] / ec ;
				pt[2] = pt[2] / ec ;
				if ( pt[0] < st[0] || pt[1] < st[1] || pt[2] < st[2] ||
					pt[0] > st[0] + len || pt[1] > st[1] + len || pt[2] > st[2] + len )
				{
					//printf("Out! %f %f %f, Box: (%d %d %d) Len: %d ec: %d\n", pt[0], pt[1], pt[2], st[0],st[1],st[2],len,ec) ;
				}

				unsigned char sg = 0 ;
				for ( int i = 0 ; i < 8 ; i ++ ) {
					if ( signs[i] == 1 ) {
						sg |= ( 1 << i ) ;
					}
					else if ( signs[i] == -1 ) {
						// Undetermined, use center sign instead
						if ( midsign == 1 ) {
							sg |= ( 1 << i ) ;
						}
						else if ( midsign == -1 ) {
							printf("Wrong!");
						}
					}
				}

				// Solve
				float mat[10];
				BoundingBoxf * box = new BoundingBoxf();
				box->begin.x = (float) st[0] ;
				box->begin.y = (float) st[1] ;
				box->begin.z = (float) st[2] ;
				box->end.x = (float) st[0] + len ;
				box->end.y = (float) st[1] + len ;
				box->end.z = (float) st[2] + len ;
				
				float error = calcPoint( ata, atb, btb, pt, mp, box, mat ) ;
#ifdef CLAMP
				if ( mp[0] < st[0] || mp[1] < st[1] || mp[2] < st[2] ||
					mp[0] > st[0] + len || mp[1] > st[1] + len || mp[2] > st[2] + len )
				{
					mp[0] = pt[0] ;
					mp[1] = pt[1] ;
					mp[2] = pt[2] ;
				}			
#endif
				if ( error <= thresh )
				{
					//mp[0] = st[0] + len / 2 ;
					//mp[1] = st[1] + len / 2 ;
					//mp[2] = st[2] + len / 2 ;

					PseudoLeafNode* pnode = new PseudoLeafNode( ht+1, sg, ata, atb, btb, mp ) ;
					for ( int i = 0 ; i < 8 ; i ++ )
						pnode->child[i] = inode->child[i] ;

					delete inode ;
					return pnode ;
				}
				else {
					return node ;
				}
			}
			
		}
		else {
			return node ;
		}
	}
	else {
		return node ;
	}
}


void Octree::readSOG( char* fname ) {
	FILE* fin = fopen( fname, "rb" ) ;
	if ( fin == NULL ) {
		printf("Can not open file %s.\n", fname) ;
	}
	
	// Process header
	char header[]="SOG.Format 1.0";
	float origin[3];
	float range;

	fread( header, sizeof( char ), strlen(header) + 1, fin );
	if ( strcmp(header, "SOG.Format 1.0") )
	{
		printf("Incorrect file format.\n", fname) ;
		exit(1);
	}
	fread( origin, sizeof( float ), 3, fin );
	fread( &range, sizeof( float ), 1, fin ) ;
	printf("Origin: (%f, %f, %f), Range: %f.\n", origin[0], origin[1], origin[2], range);

	int nlen = 128 - 4 * 4 - strlen(header) - 1 ;
	char* header2 = new char[ nlen ];
	fread( header2, sizeof( char ), nlen, fin ) ;


	fread( &(this->dimen), sizeof( int ), 1, fin ) ;
	this->maxDepth = 0 ;
	int temp = 1 ;
	while ( temp < this->dimen )
	{
		maxDepth ++ ;
		temp <<= 1 ;
	}
	printf("Dimensions: %d Depth: %d\n", this->dimen, maxDepth ) ;

	// Recursive reader
	int st[3]={0,0,0};
	this->root = readSOG( fin, st, dimen, maxDepth, origin, range ) ;
	
	printf("Done reading.\n") ;	
	fclose( fin ) ;
}

OctreeNode* Octree::readSOG( FILE* fin, int st[3], int len, int ht, float origin[3], float range )
{
	OctreeNode* rvalue = NULL ;
	int i ;

	// Get type
	char type ;
	fread( &type, sizeof( char ), 1, fin ) ;

	if ( type == 0 )
	{
		// Internal node
		rvalue = new InternalNode( ) ;

		int nlen = len / 2 ;
		int nst[3] ;

		for ( i = 0 ; i < 8 ; i ++ )
		{
			nst[0] = st[0] + vertMap[i][0] * nlen ;
			nst[1] = st[1] + vertMap[i][1] * nlen ;
			nst[2] = st[2] + vertMap[i][2] * nlen ;

			((InternalNode *)rvalue)->child[i] = readSOG( fin, nst, nlen, ht - 1, origin, range ) ;
		}
	}
	else if ( type == 1 )
	{
		// Empty node
		char sg ;
		fread( &sg, sizeof( char ), 1, fin ) ;
		sg = ~sg;
	}
	else if ( type == 2 )
	{
		// Leaf node
		unsigned char sg ;
		fread( &sg, sizeof( unsigned char ), 1, fin ) ;
		sg = ~sg;
		
		float coord[3] ;
		fread( coord, sizeof( float ), 3, fin ) ;

#ifdef SOG_RELATIVE
		coord[0] = st[0] + len * coord[0] ;
		coord[1] = st[1] + len * coord[1] ;
		coord[2] = st[2] + len * coord[2] ;
#else
		for ( i = 0 ; i < 3 ; i ++ )
		{
			coord[ i ] = ( coord[ i ] - origin[ i ] ) * dimen / range ;
		}
#endif

		rvalue = new LeafNode( ht, sg, coord ) ;
	}
	else if ( type == 3 )
	{
		// Pseudo-leaf node
		unsigned char sg ;
		fread( &sg, sizeof( unsigned char ), 1, fin ) ;
		sg = ~sg;
		
		float coord[3] ;
		fread( coord, sizeof( float ), 3, fin ) ;
		
#ifdef SOG_RELATIVE
		coord[0] = st[0] + len * coord[0] ;
		coord[1] = st[1] + len * coord[1] ;
		coord[2] = st[2] + len * coord[2] ;
#else
		for ( i = 0 ; i < 3 ; i ++ )
		{
			coord[ i ] = ( coord[ i ] - origin[ i ] ) * dimen / range ;
		}
#endif
		
		rvalue = new PseudoLeafNode( ht, sg, coord ) ;
		
		int nlen = len / 2 ;
		int nst[3] ;

		for ( int k = 0 ; i < 8 ; i ++ )
		{
			nst[0] = st[0] + vertMap[i][0] * nlen ;
			nst[1] = st[1] + vertMap[i][1] * nlen ;
			nst[2] = st[2] + vertMap[i][2] * nlen ;
			((PseudoLeafNode *)rvalue)->child[i] = readSOG( fin, nst, nlen, ht - 1, origin, range ) ;
		}
	}
	else
	{
		printf("Wrong! Type: %d\n", type);
	}
	return rvalue ;
}

void Octree::readDCF( char* fname ) {
	FILE* fin = fopen( fname, "rb" ) ;
	if ( fin == NULL )
		printf("Can not open file %s.\n", fname) ;
	
	// Process header
	char version[10] ;
	fread( version, sizeof( char ), 10, fin ) ;
	if ( strcmp( version, "multisign" ) != 0 ) {
		printf("Wrong DCF version.\n") ;
		exit(0) ;
	}
	
	fread( &(this->dimen), sizeof( int ), 1, fin ) ;
	fread( &(this->dimen), sizeof( int ), 1, fin ) ;
	fread( &(this->dimen), sizeof( int ), 1, fin ) ;
	this->maxDepth = 0 ;
	int temp = 1 ;
	while ( temp < this->dimen ) {
		maxDepth ++ ;
		temp <<= 1 ;
	}
	printf("Dimensions: %d Depth: %d\n", this->dimen, maxDepth ) ;

	// Recursive reader
	int st[3] = {0, 0, 0} ;
	this->root = readDCF( fin, st, dimen, maxDepth ) ;
	
	printf("Done reading.\n") ;	
	fclose( fin ) ;
}

OctreeNode* Octree::readDCF( FILE* fin, int st[3], int len, int ht ) {
	OctreeNode* rvalue = NULL ;

	int type ;
	fread( &type, sizeof( int ), 1, fin ) ;  // Get type
	//printf("Type: %d\n", type);

	if ( type == 0 ) { // Internal node
		rvalue = new InternalNode() ;
		int nlen = len / 2 ; // len of child node is half that of parent
		int nst[3] ;

		for ( int i = 0 ; i < 8 ; i ++ ) {  // create eight child nodes
			nst[0] = st[0] + vertMap[i][0] * nlen; // child st position
			nst[1] = st[1] + vertMap[i][1] * nlen;
			nst[2] = st[2] + vertMap[i][2] * nlen;
			((InternalNode *)rvalue)->child[i] = readDCF( fin, nst, nlen, ht - 1 ) ; // height is one less than parent
		}
	}
	else if ( type == 1 ) { // Empty node
		short sg ;
		fread( &sg, sizeof( short ), 1, fin ) ; // signs not used??
	}
	else if ( type == 2 ) { // Leaf node
		short rsg[8] ;
		fread( rsg, sizeof( short ), 8, fin ) ;
		unsigned char sg = 0 ;
		for ( int i = 0 ; i < 8 ; i ++ ) { // set signs
			if ( rsg[i] != 0 ) {
				sg |= ( 1 << i ) ;
			}
		}

		// intersections and normals
		float inters[12][3], norms[12][3] ;
		int numinters = 0;
		for ( int i = 0 ; i < 12 ; i ++ ) { // potentially there are 12 intersections
			int num ;
			fread( &num, sizeof( int ), 1, fin ) ;
			if ( num > 0 ) {
				for ( int j = 0 ; j < num ; j ++ ) {
					float off ;
					fread( &off, sizeof( float ), 1, fin ) ;
					fread( norms[numinters], sizeof( float ), 3, fin ) ;
					int dir = i / 4 ;
					int base = edgevmap[ i ][ 0 ] ;
					inters[numinters][0] = st[0] + vertMap[base][0] * len ;
					inters[numinters][1] = st[1] + vertMap[base][1] * len ;
					inters[numinters][2] = st[2] + vertMap[base][2] * len ;
					inters[numinters][dir] += off ;
					numinters ++ ;
				}
			}
		}
		
		if ( numinters > 0 ) {
			rvalue = new LeafNode( ht, (unsigned char)sg, st, len, numinters, inters, norms ) ;
		}
		else {
			rvalue = NULL ;
		}
	}
	else {
		printf("Wrong! Type: %d\n", type);
	}
	return rvalue ;
}

// no-intersections algorithm
// fname is the PLY output file
void Octree::genContourNoInter2( char* fname ) {
	int numTris = 0 ;
	int numVertices = 0 ;
	IndexedTriangleList* tlist = new IndexedTriangleList();
	VertexList* vlist = new VertexList();
	tlist->next = NULL ;
	vlist->next = NULL ;

	founds = 0 ;
	news = 0 ;

	// generate triangles
	faceVerts = 0 ;
	edgeVerts = 0 ;
	HashMap* hash = new HashMap();
	int st[3] = {0,0,0};
	printf("Processing contour...\n") ;

	clock_t start = clock( ) ;
	// one cellProc call to root processes entire tree
	cellProcContourNoInter2( root, st, dimen, hash, tlist, numTris, vlist, numVertices ) ;
	clock_t finish = clock( ) ;
	printf("Time used: %f seconds.\n", (float) (finish - start) / (float) CLOCKS_PER_SEC ) ;
	
	printf("Face vertices: %d Edge vertices: %d\n", faceVerts, edgeVerts ) ;
	printf("New hash entries: %d. Found times: %d\n", news, founds) ;

	// Finally, turn into PLY
	FILE* fout = fopen ( fname, "wb" ) ;
	printf("Vertices counted: %d Triangles counted: %d \n", numVertices, numTris ) ;
	PLYWriter::writeHeader( fout, numVertices, numTris ) ;

	VertexList* v = vlist->next ;
	while ( v != NULL ) {
		PLYWriter::writeVertex( fout, v->vt ) ;
		v = v->next ;
	}

	IndexedTriangleList* t = tlist->next ;
	for ( int i = 0 ; i < numTris ; i ++ ) {
		int inds[] = {numVertices - 1 - t->vt[0], numVertices - 1 - t->vt[1], numVertices - 1 - t->vt[2]} ;
		PLYWriter::writeFace( fout, 3, inds ) ;
		t = t->next ;
	}

	fclose( fout ) ;

	// Clear up
	delete hash ;
	v = vlist ;
	while ( v != NULL ) {
		vlist = v->next ;
		delete v ;
		v = vlist ;
	}
	t = tlist ;
	while ( t != NULL ) {
		tlist = t->next ;
		delete t ;
		t = tlist ;
	}
}


// original algorithm, may produce intersecting polygons?
void Octree::genContour( char* fname ) {
	int numTris = 0 ;
	int numVertices = 0 ;

	FILE* fout = fopen ( fname, "wb" ) ;
	cellProcCount ( root, numVertices, numTris ) ;
	printf("Vertices counted: %d Triangles counted: %d \n", numVertices, numTris ) ;
	PLYWriter::writeHeader( fout, numVertices, numTris ) ;
	int offset = 0; // ?

	clock_t start = clock( ) ;
	generateVertexIndex( root, offset, fout ) ; // write vertices to file, populate node->index
	actualTris = 0 ;
	cellProcContour( this->root, fout ) ; // a single call to root runs algorithm on entire tree
	clock_t finish = clock( ) ;
	printf("Time used: %f seconds.\n", (float) (finish - start) / (float) CLOCKS_PER_SEC ) ;
	printf("Actual triangles written: %d\n", actualTris ) ;
	fclose( fout ) ;
}

// this writes out octree vertices to the PLY file
// each vertex gets an index, which is stored in node->index
void Octree::generateVertexIndex( OctreeNode* node, int& offset, FILE* fout ) {
	NodeType type = node->getType() ;

	if ( type == INTERNAL ) { // Internal node, recurse into tree
		InternalNode* inode = ( (InternalNode* ) node ) ;
		for ( int i = 0 ; i < 8 ; i ++ ) {
			if ( inode->child[i] != NULL )
				generateVertexIndex( inode->child[i], offset, fout ) ;
		}
	}
	else if ( type == LEAF ) { // Leaf node
		LeafNode* lnode = ((LeafNode *) node) ;
		PLYWriter::writeVertex( fout, lnode->mp ) ; // write out mp (?)
		lnode->index = offset;
		offset++;
	}
	else if ( type == PSEUDOLEAF ) { // Pseudo leaf node
		PseudoLeafNode* pnode = ((PseudoLeafNode *) node) ;
		PLYWriter::writeVertex( fout, pnode->mp ) ; // write out mp
		pnode->index = offset;
		offset++;
	}
}

void Octree::cellProcContour( OctreeNode* node, FILE* fout )  {
	if ( node == NULL )
		return ;

	int type = node->getType() ;

	if ( type == 0 ) { // internal node
		InternalNode* inode = (( InternalNode * ) node ) ; //why cast?

		// 8 Cell calls on children
		for ( int i = 0 ; i < 8 ; i ++ ) 
			cellProcContour( inode->child[ i ], fout );

		// 12 face calls, faces between each child node
		for ( int i = 0 ; i < 12 ; i ++ ) {
			int c[ 2 ] = { cellProcFaceMask[ i ][ 0 ], cellProcFaceMask[ i ][ 1 ] };
			OctreeNode* fcd[2];
			fcd[0] = inode->child[ c[0] ] ;
			fcd[1] = inode->child[ c[1] ] ;
			faceProcContour( fcd, cellProcFaceMask[ i ][ 2 ], fout ) ;
		}

		// 6 edge calls
		for ( int i = 0 ; i < 6 ; i ++ ) {
			int c[ 4 ] = { cellProcEdgeMask[ i ][ 0 ], cellProcEdgeMask[ i ][ 1 ], cellProcEdgeMask[ i ][ 2 ], cellProcEdgeMask[ i ][ 3 ] };
			OctreeNode* ecd[4] ;
			for ( int j = 0 ; j < 4 ; j ++ )
				ecd[j] = inode->child[ c[j] ] ;

			edgeProcContour( ecd, cellProcEdgeMask[ i ][ 4 ], fout ) ;
		}
	}
};

// node[2] are the two nodes that share a face
// dir comes from cellProcFaceMask[i][2]  where i=0..11
void Octree::faceProcContour ( OctreeNode* node[2], int dir, FILE* fout ) 
{
	// printf("I am at a face! %d\n", dir ) ;
	if ( ! ( node[0] && node[1] ) ) {
		// printf("I am none.\n") ;
		return ;
	}

	NodeType type[2] = { node[0]->getType(), node[1]->getType() } ;

	if ( type[0] == INTERNAL || type[1] == INTERNAL ) { // both nodes internal
		// 4 face calls
		OctreeNode* fcd[2] ;
		for ( int i = 0 ; i < 4 ; i ++ ) {
			int c[2] = { faceProcFaceMask[ dir ][ i ][ 0 ], faceProcFaceMask[ dir ][ i ][ 1 ] };
			for ( int j = 0 ; j < 2 ; j ++ ) {
				if ( type[j] > 0 )
					fcd[j] = node[j];
				else 
					fcd[j] = ((InternalNode *) node[ j ] )->child[ c[j] ];
			}
			faceProcContour( fcd, faceProcFaceMask[ dir ][ i ][ 2 ], fout ) ;
		}

		// 4 edge calls
		int orders[2][4] = {{ 0, 0, 1, 1 }, { 0, 1, 0, 1 }} ;
		OctreeNode* ecd[4] ;
			
		for ( int i = 0 ; i < 4 ; i ++ ) {
			int c[4] = { faceProcEdgeMask[ dir ][ i ][ 1 ], faceProcEdgeMask[ dir ][ i ][ 2 ],
						 faceProcEdgeMask[ dir ][ i ][ 3 ], faceProcEdgeMask[ dir ][ i ][ 4 ] };
			int* order = orders[ faceProcEdgeMask[ dir ][ i ][ 0 ] ] ;

			for ( int j = 0 ; j < 4 ; j ++ ) {
				if ( type[order[j]] > 0 )
					ecd[j] = node[order[j]] ;
				else
					ecd[j] = ( (InternalNode *) node[ order[ j ] ] )->child[ c[j] ] ;
			}
			edgeProcContour( ecd, faceProcEdgeMask[ dir ][ i ][ 5 ], fout ) ;
		}
//		printf("I am done.\n") ;
	}
	else {
//		printf("I don't have any children.\n") ;
	}
};

// a common edge between four nodes in node[4]
// "dir" comes from cellProcEdgeMask
void Octree::edgeProcContour ( OctreeNode* node[4], int dir, FILE* fout ) 
{
	if ( ! ( node[0] && node[1] && node[2] && node[3] ) )
		return ;


	NodeType type[4] = { node[0]->getType(), node[1]->getType(), node[2]->getType(), node[3]->getType() } ;

	if ( type[0] != INTERNAL && type[1] != INTERNAL  && type[2] != INTERNAL && type[3] != INTERNAL ) {
		processEdgeWrite( node, dir, fout ) ; // a quad/traingle is output?
	}
	else {
		// 2 edge calls
		OctreeNode* ecd[4] ;
		for ( int i = 0 ; i < 2 ; i ++ ) {
			int c[ 4 ] = { edgeProcEdgeMask[ dir ][ i ][ 0 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 1 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 2 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 3 ] } ;

			for ( int j = 0 ; j < 4 ; j ++ ) {
				if ( type[j] > 0 )
					ecd[j] = node[j] ;
				else
					ecd[j] = ((InternalNode *) node[j])->child[ c[j] ] ;
			}

			edgeProcContour( ecd, edgeProcEdgeMask[ dir ][ i ][ 4 ], fout ) ;
		}

	}
};

void Octree::processEdgeWrite ( OctreeNode* node[4], int dir, FILE* fout )  {
	// Get minimal cell
	int type, ht, minht = this->maxDepth+1, mini = -1 ;
	int ind[4], sc[4], flip[4] = {0,0,0,0} ;
	int flip2;
	for ( int i = 0 ; i < 4 ; i ++ ) {
		if ( node[i]->getType() == LEAF ) {
			LeafNode* lnode = ((LeafNode *) node[i]) ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( lnode->height < minht ) {
				minht = lnode->height ;
				mini = i ;
				if ( lnode->getSign(c1) > 0 )
				{
					flip2 = 1 ;
				}
				else
				{
					flip2 = 0 ;
				}
			}
			ind[i] = lnode->index ;

			if ( lnode->getSign( c1 ) == lnode->getSign( c2 ) )
				sc[ i ] = 0 ;
			else
				sc[ i ] = 1 ;

//				if ( lnode->getSign(c1) > 0 )
//				{
//					flip[ i ] = 1 ;
//				}
			// }
		}
		else {
			PseudoLeafNode* pnode = ((PseudoLeafNode *) node[i]) ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( pnode->height < minht )
			{
				minht = pnode->height ;
				mini = i ;
				if ( pnode->getSign(c1) > 0 )
					flip2 = 1 ;
				else
					flip2 = 0 ;
			}
			ind[i] = pnode->index ;

			if ( pnode->getSign( c1 ) == pnode->getSign( c2 ) )
				sc[ i ] = 0 ;
			else
				sc[ i ] = 1 ;

//				if ( pnode->getSign(c1) > 0 )
//				{
//					flip[ i ] = 1 ;
//					flip2 = 1 ;
//				}
//			}
		}

	}

	if ( sc[ mini ] == 1 )
	{
		if ( flip2 == 0 )
		{
			actualTris ++ ;
			if ( ind[0] == ind[1] )
			{
				int tind[] = { ind[0], ind[3], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else if ( ind[1] == ind[3] )
			{
				int tind[] = { ind[0], ind[1], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else if ( ind[3] == ind[2] )
			{
				int tind[] = { ind[0], ind[1], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else if ( ind[2] == ind[0] )
			{
				int tind[] = { ind[1], ind[3], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else
			{
				int tind1[] = { ind[0], ind[1], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind1 ) ;
				int tind2[] = { ind[0], ind[3], ind[2] } ;
				PLYWriter::writeFace( fout, 3, tind2 ) ;
				actualTris ++ ;
			}
		}
		else
		{
			actualTris ++ ;
			if ( ind[0] == ind[1] )
			{
				int tind[] = { ind[0], ind[2], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else if ( ind[1] == ind[3] )
			{
				int tind[] = { ind[0], ind[2], ind[1] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else if ( ind[3] == ind[2] )
			{
				int tind[] = { ind[0], ind[3], ind[1] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else if ( ind[2] == ind[0] )
			{
				int tind[] = { ind[1], ind[2], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind ) ;
			}
			else
			{
				int tind1[] = { ind[0], ind[3], ind[1] } ;
				PLYWriter::writeFace( fout, 3, tind1 ) ;
				int tind2[] = { ind[0], ind[2], ind[3] } ;
				PLYWriter::writeFace( fout, 3, tind2 ) ;
				actualTris ++ ;
			}
		}

	}

};


void Octree::cellProcCount( OctreeNode* node, int& nverts, int& nfaces )  {
	if ( node == NULL )
		return ;

	int type = node->getType() ;

	if ( type > 0 ) {
		nverts ++ ;
	}
	else
	{
		InternalNode* inode = (( InternalNode * ) node ) ;

		// 8 Cell calls
		for ( int i = 0 ; i < 8 ; i ++ )
			cellProcCount( inode->child[ i ], nverts, nfaces ) ;

		// 12 face calls
		OctreeNode* fcd[2] ;
		for ( int i = 0 ; i < 12 ; i ++ ) {
			int c[ 2 ] = { cellProcFaceMask[ i ][ 0 ], cellProcFaceMask[ i ][ 1 ] };

			fcd[0] = inode->child[ c[0] ] ;
			fcd[1] = inode->child[ c[1] ] ;

			faceProcCount( fcd, cellProcFaceMask[ i ][ 2 ], nverts, nfaces ) ;
		}

		// 6 edge calls
		OctreeNode* ecd[4] ;
		for ( int i = 0 ; i < 6 ; i ++ )
		{
			int c[ 4 ] = { cellProcEdgeMask[ i ][ 0 ], cellProcEdgeMask[ i ][ 1 ], cellProcEdgeMask[ i ][ 2 ], cellProcEdgeMask[ i ][ 3 ] };

			for ( int j = 0 ; j < 4 ; j ++ )
				ecd[j] = inode->child[ c[j] ] ;

			edgeProcCount( ecd, cellProcEdgeMask[ i ][ 4 ], nverts, nfaces ) ;
		}
	}
};

void Octree::faceProcCount ( OctreeNode* node[2], int dir, int& nverts, int& nfaces ) 
{
	if ( ! ( node[0] && node[1] ) )
		return ;

	int type[2] = { node[0]->getType(), node[1]->getType() } ;

	if ( type[0] == 0 || type[1] == 0 )
	{
		int i, j ;

		// 4 face calls
		OctreeNode* fcd[2] ;
		for ( int i = 0 ; i < 4 ; i ++ )
		{
			int c[2] = { faceProcFaceMask[ dir ][ i ][ 0 ], faceProcFaceMask[ dir ][ i ][ 1 ] };
			for ( int j = 0 ; j < 2 ; j ++ )
			{
				if ( type[j] > 0 )
				{
					fcd[j] = node[j] ;
				}
				else
				{
					fcd[j] = ((InternalNode *) node[ j ] )->child[ c[j] ] ;
				}
			}
			faceProcCount( fcd, faceProcFaceMask[ dir ][ i ][ 2 ], nverts, nfaces ) ;
		}

		// 4 edge calls
		int orders[2][4] = {{ 0, 0, 1, 1 }, { 0, 1, 0, 1 }} ;
		OctreeNode* ecd[4] ;
			
		for ( int i = 0 ; i < 4 ; i ++ )
		{
			int c[4] = { faceProcEdgeMask[ dir ][ i ][ 1 ], faceProcEdgeMask[ dir ][ i ][ 2 ],
						 faceProcEdgeMask[ dir ][ i ][ 3 ], faceProcEdgeMask[ dir ][ i ][ 4 ] };
			int* order = orders[ faceProcEdgeMask[ dir ][ i ][ 0 ] ] ;

			for ( int j = 0 ; j < 4 ; j ++ )
			{
				if ( type[order[j]] > 0 )
				{
					ecd[j] = node[order[j]] ;
				}
				else
				{
					ecd[j] = ( (InternalNode *) node[ order[ j ] ] )->child[ c[j] ] ;
				}
			}

			edgeProcCount( ecd, faceProcEdgeMask[ dir ][ i ][ 5 ], nverts, nfaces ) ;
		}
	}
};

void Octree::edgeProcCount ( OctreeNode* node[4], int dir, int& nverts, int& nfaces ) 
{
	if ( ! ( node[0] && node[1] && node[2] && node[3] ) )
	{
		return ;
	}

	int type[4] = { node[0]->getType(), node[1]->getType(), node[2]->getType(), node[3]->getType() } ;

	if ( type[0] > 0 && type[1] > 0 && type[2] > 0 && type[3] > 0 )
	{
		processEdgeCount( node, dir, nverts, nfaces ) ;
	}
	else
	{
		int i, j ;

		// 2 edge calls
		OctreeNode* ecd[4] ;
		for ( int i = 0 ; i < 2 ; i ++ )
		{
			int c[ 4 ] = { edgeProcEdgeMask[ dir ][ i ][ 0 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 1 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 2 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 3 ] } ;

			for ( int j = 0 ; j < 4 ; j ++ )
			{
				if ( type[j] > 0 )
				{
					ecd[j] = node[j] ;
				}
				else
				{
					ecd[j] = ((InternalNode *) node[j])->child[ c[j] ] ;
				}
			}

			edgeProcCount( ecd, edgeProcEdgeMask[ dir ][ i ][ 4 ], nverts, nfaces ) ;
		}

	}
};

void Octree::processEdgeCount ( OctreeNode* node[4], int dir, int& nverts, int& nfaces ) 
{
	// Get minimal cell
	int i, type, ht, minht = maxDepth+1, mini = -1 ;
	int ind[4], sc[4], flip[4] = {0,0,0,0} ;
	for ( i = 0 ; i < 4 ; i ++ ) {
		if ( node[i]->getType() == 1 ) {
			LeafNode* lnode = ((LeafNode *) node[i]) ;

			if ( lnode->height < minht ) {
				minht = lnode->height ;
				mini = i ;
			}
			ind[i] = lnode->index ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( lnode->getSign( c1 ) == lnode->getSign( c2 ) ) {
				sc[ i ] = 0 ;
			}
			else {
				sc[ i ] = 1 ;
				if ( lnode->getSign(c1) > 0 )
					flip[ i ] = 1 ;
			}
		}
		else {
			PseudoLeafNode* pnode = ((PseudoLeafNode *) node[i]) ;
			if ( pnode->height < minht ) {
				minht = pnode->height ;
				mini = i ;
			}
			ind[i] = pnode->index ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( pnode->getSign( c1 ) == pnode->getSign( c2 ) )
			{
				sc[ i ] = 0 ;
			}
			else
			{
				sc[ i ] = 1 ;

				if ( pnode->getSign(c1) > 0 )
				{
					flip[ i ] = 1 ;
				}
			}
		}

	}

	if ( sc[ mini ] == 1 ) {
		nfaces ++ ;
		if ( node[0] != node[1] && node[1] != node[3] && node[3] != node[2] && node[2] != node[0] )
		{
			nfaces ++ ;
		}

	}

};

/************************************************************************/
/* Start Non-inters                                                     */
/************************************************************************/



int Octree::testFace( int st[3], int len, int dir, float v1[3], float v2[3] ) 
{
#ifdef TESS_UNIFORM
	return 0 ;
#endif

#ifdef TESS_NONE
	return 1 ;
#endif
	float vec[3] = { v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2] } ;
	float ax1[3], ax2[3] ;
	float ed1[3]={0,0,0}, ed2[3] = {0,0,0};
	ed1[(dir+1)%3]=1;
	ed2[(dir+2)%3]=1;

	Intersection::cross( ed1, vec, ax1 ) ;
	Intersection::cross( ed2, vec, ax2 ) ;

	Triangle* t1 = new Triangle ;
	Triangle* t2 = new Triangle ;

	for ( int i = 0 ; i < 3 ; i ++ )
	{
		t1->vt[0][i] = v1[i] ;
		t1->vt[1][i] = v2[i] ;
		t1->vt[2][i] = v2[i] ;

		t2->vt[0][i] = st[i] ;
		t2->vt[1][i] = st[i] ;
		t2->vt[2][i] = st[i] ;
	}
	t2->vt[1][(dir+1)%3] += len ;
	t2->vt[2][(dir+2)%3] += len ;

	if ( Intersection::separating( ax1, t1, t2 ) || Intersection::separating( ax2, t1, t2 ) )
	{
		faceVerts ++ ;
/*
			printf("\n{{{%d, %d, %d},%d,%d},{{%f, %f, %f},{%f, %f, %f}}}\n",
				st[0],st[1], st[2],
				len, dir+1,
				v1[0], v1[1], v1[2],
				v2[0], v2[1], v2[2]) ;
*/		
		return 0 ;
	}
	else
	{
		return 1 ;
	}
};

int Octree::testEdge( int st[3], int len, int dir, OctreeNode* node[4], float v[4][3] ) 
{
#ifdef TESS_UNIFORM
	return 0 ;
#endif

#ifdef TESS_NONE
	return 1 ;
#endif

	if ( node[0] == node[1] || node[1] == node[3] || node[3] == node[2] || node[2] == node[0] )
	{
//		return 1 ;
	}

	float p1[3] ={st[0], st[1], st[2]};
	float p2[3] ={st[0], st[1], st[2]};
	p2[dir] += len ;

	int nbr[]={0,1,3,2,0} ;
	int nbr2[]={3,2,0,1,3} ;
	float nm[3], vec1[4][3], vec2[4][3] ;
	float d1, d2 ;

	for ( int i = 0 ; i < 4 ; i ++ )
	{
		for ( int j = 0 ; j < 3 ; j ++ )
		{
			vec1[i][j] = v[i][j] - p1[j] ;
			vec2[i][j] = v[i][j] - p2[j] ;
		}
	}

#ifdef EDGE_TEST_CONVEXITY
	for ( i = 0 ; i < 4 ; i ++ )
	{
		int a = nbr[i] ;
		int b = nbr[i+1] ;
		int c = nbr2[i] ;
		int d = nbr2[i+1] ;

		if ( node[a] == node[b] )
		{
			continue ;
		}

		Intersection::cross( vec1[a], vec1[b], nm ) ;
		d1 = Intersection::dot( vec1[c], nm ) ;
		d2 = Intersection::dot( vec1[d], nm ) ;

		if ( d1 * d2 < 0 )
		{
			/*
			printf("1\n{{{%f, %f, %f},{%f, %f, %f}},{{%f, %f, %f},{%f, %f, %f},{%f, %f, %f},{%f, %f, %f}}}\n",
				p1[0], p1[1], p1[2],
				p2[0], p2[1], p2[2],
				v[0][0], v[0][1], v[0][2],
				v[1][0], v[1][1], v[1][2],
				v[2][0], v[2][1], v[2][2],
				v[3][0], v[3][1], v[3][2]) ;
			*/
			return 0 ;
		}

		Intersection::cross( vec2[a], vec2[b], nm ) ;
		d1 = Intersection::dot( vec2[c], nm ) ;
		d2 = Intersection::dot( vec2[d], nm ) ;

		if ( d1 * d2 < 0 )
		{
			/*
			printf("2\n{{{%f, %f, %f},{%f, %f, %f}},{{%f, %f, %f},{%f, %f, %f},{%f, %f, %f},{%f, %f, %f}}}\n",
				p1[0], p1[1], p1[2],
				p2[0], p2[1], p2[2],
				v[0][0], v[0][1], v[0][2],
				v[1][0], v[1][1], v[1][2],
				v[2][0], v[2][1], v[2][2],
				v[3][0], v[3][1], v[3][2]) ;
			*/
			return 0 ;
		}
	}
#else
#ifdef EDGE_TEST_FLIPDIAGONAL
	Triangle* t1 = new Triangle ;
	Triangle* t2 = new Triangle ;
	int tri[2][2][4] = {{{0,1,3,2},{3,2,0,1}},{{2,0,1,3},{1,3,2,0}}} ;
	for ( i = 0 ; i < 2 ; i ++ )
	{
		int good = 1 ;
		for ( int j = 0 ; j < 2 ; j ++ )
		{
			// Top
			for ( int k = 0 ; k < 3 ; k ++ )
			{
				t1->vt[0][k] = v[tri[i][j][0]][k] ;
				t1->vt[1][k] = v[tri[i][j][1]][k] ;
				t1->vt[2][k] = v[tri[i][j][2]][k] ;

				t2->vt[0][k] = p1[k] ;
				t2->vt[1][k] = v[tri[i][j][2]][k] ;
				t2->vt[2][k] = v[tri[i][j][3]][k] ;
			}

			if ( Intersection::testIntersection( t1, t2 ) )
			{
				good = 0 ;
				break ;
			}
			
			// Bottom
			for ( k = 0 ; k < 3 ; k ++ )
			{
				t2->vt[0][k] = p2[k] ;
				t2->vt[1][k] = v[tri[i][j][2]][k] ;
				t2->vt[2][k] = v[tri[i][j][3]][k] ;
			}
			
			if ( Intersection::testIntersection( t1, t2 ) )
			{
				good = 0 ;
				break ;
			}
		}

		if ( good )
		{
			return (i+1) ;
		}
	}
	return 0 ;

#else
#ifdef EDGE_TEST_NEW
	Triangle* t1 = new Triangle ;
	int tri[2][2][4] = {{{0,1,3,2},{3,2,0,1}},{{2,0,1,3},{1,3,2,0}}} ;
	for ( int i = 0 ; i < 2 ; i ++ )
	{
		// For each triangulation
		for ( int j = 0 ; j < 2 ; j ++ )
		{
			// Starting with each triangle
			int k ;

			// Check triangle and dual edge
			float vec1[3], vec2[3], vec3[3], vec[3] ;
			for ( k = 0 ; k < 3 ; k ++ )
			{
				t1->vt[0][k] = v[tri[i][j][0]][k] ;
				t1->vt[1][k] = v[tri[i][j][1]][k] ;
				t1->vt[2][k] = v[tri[i][j][2]][k] ;

				vec1[k] = t1->vt[1][k] - t1->vt[0][k] ;
				vec2[k] = t1->vt[2][k] - t1->vt[1][k] ;
			}

			float axes[3] ;
			Intersection::cross( vec1, vec2, axes ) ;

			if ( Intersection::separating( axes, t1, p1, p2 ) )
			{
				continue ;
			}
			
			// Check diagonal and the other triangle
			for ( k = 0 ; k < 3 ; k ++ )
			{
				t1->vt[0][k] = p1[k] ;
				t1->vt[1][k] = p2[k] ;
				t1->vt[2][k] = v[tri[i][j][3]][k] ;

				vec1[k] = t1->vt[1][k] - t1->vt[0][k] ;
				vec2[k] = t1->vt[2][k] - t1->vt[1][k] ;
				vec3[k] = t1->vt[0][k] - t1->vt[2][k] ;

				vec[k] = v[tri[i][j][2]][k] - v[tri[i][j][0]][k] ;
			}

			float axes1[3], axes2[3], axes3[3] ;
			Intersection::cross( vec1, vec, axes1 ) ;
			Intersection::cross( vec2, vec, axes2 ) ;
			Intersection::cross( vec3, vec, axes3 ) ;

			if ( Intersection::separating( axes1, t1, v[tri[i][j][0]], v[tri[i][j][2]] ) || 
				 Intersection::separating( axes2, t1, v[tri[i][j][0]], v[tri[i][j][2]] ) ||
				 Intersection::separating( axes3, t1, v[tri[i][j][0]], v[tri[i][j][2]] ) )
			{
				continue ;
			}

			return (i+1) ;
		}
	}
	return 0 ;


#endif

#endif
#endif
	
	return 1 ;
};


void Octree::makeEdgeVertex( int st[3], int len, int dir, OctreeNode* node[4], float mp[4][3], float v[3] ) 
{
	int nlen = len / 2 ;
	v[0] = st[0] ;
	v[1] = st[1] ;
	v[2] = st[2] ;
	v[dir] += nlen ;

	//return ;

//	if ( this->hasQEF == 0 )
//	{
//		return ;
//	}


	/* QEF based method 

	int i, j ;
	float ata[6] = { 0, 0, 0, 0, 0, 0 };
	float atb[3] = { 0, 0, 0 } ;
	float btb = 0 ;

	// Gather QEF
	for ( j = 0 ; j < 4 ; j ++ )
	{
		if ( node[j]->getType() == 1 )
		{
			LeafNode* lnode = ((LeafNode *) node[j]) ;
			for ( i = 0 ; i < 6 ; i ++ )
			{
				ata[i] += lnode->ata[i] ;
			}
			for ( i = 0 ; i < 3 ; i ++ )
			{
				atb[i] += lnode->atb[i] ;
			}
			btb += lnode->btb ;
		}
		else
		{
			PseudoLeafNode* lnode = ((PseudoLeafNode *) node[j]) ;
			for ( i = 0 ; i < 6 ; i ++ )
			{
				ata[i] += lnode->ata[i] ;
			}
			for ( i = 0 ; i < 3 ; i ++ )
			{
				atb[i] += lnode->atb[i] ;
			}
			btb += lnode->btb ;
		}
	}

	// Find coefficient A and b
	float A, b ;
	switch ( dir )
	{
	case 0 : 
		A = ata[ 0 ] ;
		b = atb[ 0 ] - ata[ 1 ] * st[ 1 ] - ata[ 2 ] * st[ 2 ] ;
		break ;
	case 1 :
		A = ata[ 3 ] ;
		b = atb[ 1 ] - ata[ 1 ] * st[ 0 ] - ata[ 4 ] * st[ 2 ] ;
		break ;
	case 2 :
		A = ata[ 5 ] ;
		b = atb[ 2 ] - ata[ 2 ] * st[ 0 ] - ata[ 4 ] * st[ 1 ] ;
		break ;
	}

	if ( A == 0 )
	{
		printf("A is zero!\n") ;
		return ;
	}
	else
	{
		v[ dir ] = b / A ;
		if ( v[ dir ] < st[ dir ] )
		{
			v[ dir ] = st[ dir ] ;
		}
		else if ( v[ dir ] > st[ dir ] + len)
		{
			v[ dir ] = st[ dir ] + len ;
		}
	}
	*/


	/* Barycentric approach */
	float pt[4][2], pv[4], x[2] ;
	int dir2 = ( dir + 1 ) % 3 ;
	int dir3 = ( dir2 + 1 ) % 3 ;
	int seq[] = {0,1,3,2} ;
	float epsilon = 0.000001f;

	for ( int i = 0 ; i < 4 ; i ++ )
	{
		pt[ i ][ 0 ] = mp[ seq[i] ][ dir2 ] ;
		pt[ i ][ 1 ] = mp[ seq[i] ][ dir3 ] ;
		pv[ i ] = mp[ seq[i] ][ dir ] ;
	}
	x[ 0 ] = st[ dir2 ] ;
	x[ 1 ] = st[ dir3 ] ;

	// Compute mean-value interpolation

	float vec[4][2], d[4] ;
	for ( int i = 0 ; i < 4 ; i ++ )
	{
		vec[i][0] = pt[i][0] - x[0] ;
		vec[i][1] = pt[i][1] - x[1] ;

		d[i] = sqrt( vec[i][0] * vec[i][0] + vec[i][1] * vec[i][1] );

		if ( d[i] < epsilon )
		{
			
			v[ dir ] = pv[ i ] ;
			if ( v[dir] < st[dir] )
			{
				v[dir] = st[dir] ;
			}
			else if ( v[dir] > st[dir] + len )
			{
				v[dir] = st[dir] + len ;
			}
			
			return ;
		}
	}

	float w[4]={0,0,0,0}, totw = 0 ;
	for ( int i = 0 ; i < 4 ; i ++ )
	{
		int i2 = ( i + 1 ) % 4 ;
		float sine = ( vec[i][0] * vec[i2][1] - vec[i][1] * vec[i2][0] ) / ( d[i] * d[i2] ) ;
		float cosine = ( vec[i][0] * vec[i2][0] + vec[i][1] * vec[i2][1] ) / ( d[i] * d[i2] ) ;

		if ( fabs(cosine + 1) < epsilon )
		{
			v[ dir ] = ( pv[ i ] * d[i2] + pv[ i2 ] * d[i] ) / (d[i] + d[i2]);
			
			if ( v[dir] < st[dir] )
			{
				v[dir] = st[dir] ;
			}
			else if ( v[dir] > st[dir] + len )
			{
				v[dir] = st[dir] + len ;
			}
			
			return ;
		}

		float tan2 = sine / ( 1 + cosine ) ;

		w[i] += ( tan2 / d[i] ) ; 
		w[i2] += ( tan2 / d[i2] ) ; 

		totw += ( tan2 / d[i] ) ;
		totw += ( tan2 / d[i2] ) ;
	}

	v[dir] = 0 ;
	for ( int i = 0 ; i < 4 ; i ++ )
	{
		v[dir] += ( w[i] * pv[ i ] / totw );
	}
	/**/
	if ( v[dir] < st[dir] )
	{
		v[dir] = st[dir] ;
	}
	else if ( v[dir] > st[dir] + len )
	{
		v[dir] = st[dir] + len ;
	}
	
};


void Octree::cellProcContourNoInter2( OctreeNode* node, int st[3], int len, 
HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts )
{
//	printf("I am at a cell! %d %d %d, %d \n", st[0], st[1], st[2], len ) ;	
	if ( node == NULL )
	{
//		printf("Empty cell.\n") ;
		return ;
	}

	int type = node->getType() ;

	if ( type == 0 )
	{
		InternalNode* inode = (( InternalNode * ) node ) ;

		// 8 Cell calls
//			printf("Process cell calls!\n") ;
		int nlen = len / 2 ;
		int nst[3] ;

		for ( int i = 0 ; i < 8 ; i ++ )
		{
//			printf("Cell %d..\n", i) ;
			nst[0] = st[0] + vertMap[i][0] * nlen ;
			nst[1] = st[1] + vertMap[i][1] * nlen ;
			nst[2] = st[2] + vertMap[i][2] * nlen ;
			cellProcContourNoInter2( inode->child[ i ], nst, nlen, hash, tlist, numTris, vlist, numVerts ) ;
//		printf("Return from %d %d %d, %d \n", nst[0], nst[1], nst[2], nlen ) ;	
		}
//					printf("I am done with cells!\n") ;

		// 12 face calls
//			printf("Process face calls!\n") ;
		OctreeNode* fcd[2] ;
		int dirCell2[3][4][3] = {
			{{0,-1,-1},{0,-1,0},{0,0,-1},{0,0,0}},
			{{-1,0,-1},{0,0,-1},{-1,0,0},{0,0,0}},
			{{-1,-1,0},{-1,0,0},{0,-1,0},{0,0,0}}};
		for ( int i = 0 ; i < 3 ; i ++ )
			for ( int j = 0 ; j < 4 ; j ++ )
			{
				nst[0] = st[0] + nlen + dirCell2[i][j][0] * nlen ;
				nst[1] = st[1] + nlen + dirCell2[i][j][1] * nlen ;
				nst[2] = st[2] + nlen + dirCell2[i][j][2] * nlen ;
				
				int ed = i * 4 + j ;
				int c[ 2 ] = { cellProcFaceMask[ ed ][ 0 ], cellProcFaceMask[ ed ][ 1 ] };
				
				fcd[0] = inode->child[ c[0] ] ;
				fcd[1] = inode->child[ c[1] ] ;
				
				faceProcContourNoInter2( fcd, nst, nlen, cellProcFaceMask[ ed ][ 2 ], hash, tlist, numTris, vlist, numVerts ) ;
			}
//					printf("I am done with faces!\n") ;

		// 6 edge calls
//			printf("Process edge calls!\n") ;
		OctreeNode* ecd[4] ;
		for ( int i = 0 ; i < 6 ; i ++ )
		{
			int c[ 4 ] = { cellProcEdgeMask[ i ][ 0 ], cellProcEdgeMask[ i ][ 1 ], cellProcEdgeMask[ i ][ 2 ], cellProcEdgeMask[ i ][ 3 ] };

			for ( int j = 0 ; j < 4 ; j ++ )
			{
				ecd[j] = inode->child[ c[j] ] ;
			}

			int dir = cellProcEdgeMask[ i ][ 4 ] ;
			nst[0] = st[0] + nlen ;
			nst[1] = st[1] + nlen ;
			nst[2] = st[2] + nlen ;
			if ( i % 2 == 0 )
			{
				nst[ dir ] -= nlen ;
			}

			edgeProcContourNoInter2( ecd, nst, nlen, dir, hash, tlist, numTris, vlist, numVerts ) ;
		}
//					printf("I am done with edges!\n") ;

	}
//	printf("I am done with cell %d %d %d, %d \n", st[0], st[1], st[2], len ) ;	
};

void Octree::faceProcContourNoInter2( OctreeNode* node[2], int st[3], int len, int dir, HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts )
{
//	printf("I am at a face! %d %d %d, %d, %d\n", st[0], st[1], st[2], len, dir ) ;
	if ( ! ( node[0] && node[1] ) )
	{
//		printf("I am none.\n") ;
		return ;
	}

	int type[2] = { node[0]->getType(), node[1]->getType() } ;

	if ( type[0] == 0 || type[1] == 0 )
	{
		int i, j ;
		int nlen = len / 2 ;
		int nst[3] ;

		// 4 face calls
		OctreeNode* fcd[2] ;
		int iface = faceProcFaceMask[ dir ][ 0 ][ 0 ] ;
		for ( i = 0 ; i < 4 ; i ++ )
		{
			int c[2] = { faceProcFaceMask[ dir ][ i ][ 0 ], faceProcFaceMask[ dir ][ i ][ 1 ] };
			for ( int j = 0 ; j < 2 ; j ++ )
			{
				if ( type[j] > 0 )
				{
					fcd[j] = node[j] ;
				}
				else
				{
					fcd[j] = ((InternalNode *) node[ j ] )->child[ c[j] ] ;
				}
			}

			nst[0] = st[0] + nlen * ( vertMap[ c[ 0 ] ][ 0 ] - vertMap[ iface ][ 0 ] );
			nst[1] = st[1] + nlen * ( vertMap[ c[ 0 ] ][ 1 ] - vertMap[ iface ][ 1 ] );
			nst[2] = st[2] + nlen * ( vertMap[ c[ 0 ] ][ 2 ] - vertMap[ iface ][ 2 ] );

			faceProcContourNoInter2( fcd, nst, nlen, faceProcFaceMask[ dir ][ i ][ 2 ], hash, tlist, numTris, vlist, numVerts ) ;
		}


		// 4 edge calls
		int orders[2][4] = {{ 0, 0, 1, 1 }, { 0, 1, 0, 1 }} ;
		OctreeNode* ecd[4] ;
			
		for ( i = 0 ; i < 4 ; i ++ )
		{
			int c[4] = { faceProcEdgeMask[ dir ][ i ][ 1 ], faceProcEdgeMask[ dir ][ i ][ 2 ],
						 faceProcEdgeMask[ dir ][ i ][ 3 ], faceProcEdgeMask[ dir ][ i ][ 4 ] };
			int* order = orders[ faceProcEdgeMask[ dir ][ i ][ 0 ] ] ;

			for ( int j = 0 ; j < 4 ; j ++ )
			{
				if ( type[order[j]] > 0 )
				{
					ecd[j] = node[order[j]] ;
				}
				else
				{
					ecd[j] = ( (InternalNode *) node[ order[ j ] ] )->child[ c[j] ] ;
				}
			}

			int ndir = faceProcEdgeMask[ dir ][ i ][ 5 ] ;
			nst[0] = st[0] + nlen ;
			nst[1] = st[1] + nlen ;
			nst[2] = st[2] + nlen ;
			nst[dir] -= nlen ;
			if ( i % 2 == 0 )
			{
				nst[ ndir ] -= nlen ;
			}

			edgeProcContourNoInter2( ecd, nst, nlen, ndir, hash, tlist, numTris, vlist, numVerts ) ;
		}
//		printf("I am done.\n") ;
	}
	else
	{
//		printf("i don't have children.\n");
	}
};

void Octree::edgeProcContourNoInter2( OctreeNode* node[4], int st[3], int len, int dir, HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts )
{
//	printf("I am at an edge! %d %d %d \n", st[0], st[1], st[2] ) ;
	if ( ! ( node[0] && node[1] && node[2] && node[3] ) )
	{
//		printf("I am done!\n") ;
		return ;
	}

	int type[4] = { node[0]->getType(), node[1]->getType(), node[2]->getType(), node[3]->getType() } ;

	if ( type[0] > 0 && type[1] > 0 && type[2] > 0 && type[3] > 0 )
	{
		this->processEdgeNoInter2( node, st, len, dir, hash, tlist, numTris, vlist, numVerts ) ;
	}
	else
	{
		int i, j ;
		int nlen = len / 2 ;
		int nst[3] ;

		// 2 edge calls
		OctreeNode* ecd[4] ;
		for ( i = 0 ; i < 2 ; i ++ )
		{
			int c[ 4 ] = { edgeProcEdgeMask[ dir ][ i ][ 0 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 1 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 2 ], 
						   edgeProcEdgeMask[ dir ][ i ][ 3 ] } ;

			for ( int j = 0 ; j < 4 ; j ++ )
			{
				if ( type[j] > 0 )
				{
					ecd[j] = node[j] ;
				}
				else
				{
					ecd[j] = ((InternalNode *) node[j])->child[ c[j] ] ;
				}
			}

			nst[0] = st[0] ;
			nst[1] = st[1] ;
			nst[2] = st[2] ;
			nst[dir] += nlen * i ;

			edgeProcContourNoInter2( ecd, nst, nlen, edgeProcEdgeMask[ dir ][ i ][ 4 ], hash, tlist, numTris, vlist, numVerts ) ;
		}

	}
//		printf("I am done!\n") ;
};


void Octree::processEdgeNoInter2( OctreeNode* node[4], int st[3], int len, int dir, HashMap* hash, IndexedTriangleList* tlist, int& numTris, VertexList* vlist, int& numVerts )
{
//	printf("I am at a leaf edge! %d %d %d\n", st[0], st[1], st[2] ) ;
	// Get minimal cell
	int i, type, minht = maxDepth+1, mini = -1 ;
	int ind[4], sc[4], ht[4], flip=0;
	float mp[4][3] ;
	for ( i = 0 ; i < 4 ; i ++ )
	{
		if ( node[i]->getType() == 1 )
		{
			LeafNode* lnode = ((LeafNode *) node[i]) ;
			ht[i] = lnode->height ;
			mp[i][0] = lnode->mp[0] ;
			mp[i][1] = lnode->mp[1] ;
			mp[i][2] = lnode->mp[2] ;


			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( lnode->height < minht )
			{
				minht = lnode->height ;
				mini = i ;
				if ( lnode->getSign(c1) > 0 )
				{
					flip = 1 ;
				}
				else
				{
					flip = 0 ;
				}
			}
			ind[i] = lnode->index ;
			if ( ind[i] < 0 )
			{
				// Create new index
				VertexList* nv = new VertexList ;
				nv->vt[0] = lnode->mp[0] ;
				nv->vt[1] = lnode->mp[1] ;
				nv->vt[2] = lnode->mp[2] ;
				nv->next = vlist->next ;
				vlist->next = nv ;
				ind[i] = numVerts ;
				lnode->index = numVerts ;
				numVerts ++ ;
			}

			if ( lnode->getSign( c1 ) == lnode->getSign( c2 ) )
			{
				sc[ i ] = 0 ;
			}
			else
			{
				sc[ i ] = 1 ;
			}
		}
		else if ( node[i]->getType() == 2 )
		{
			PseudoLeafNode* pnode = ((PseudoLeafNode *) node[i]) ;
			ht[i] = pnode->height ;
			mp[i][0] = pnode->mp[0] ;
			mp[i][1] = pnode->mp[1] ;
			mp[i][2] = pnode->mp[2] ;

			int ed = processEdgeMask[dir][i] ;
			int c1 = edgevmap[ed][0] ;
			int c2 = edgevmap[ed][1] ;

			if ( pnode->height < minht )
			{
				minht = pnode->height ;
				mini = i ;
				if ( pnode->getSign(c1) > 0 )
				{
					flip = 1 ;
				}
				else
				{
					flip = 0 ;
				}
			}
			ind[i] = pnode->index ;
			if ( ind[i] < 0 )
			{
				// Create new index
				VertexList* nv = new VertexList ;
				nv->vt[0] = pnode->mp[0] ;
				nv->vt[1] = pnode->mp[1] ;
				nv->vt[2] = pnode->mp[2] ;
				nv->next = vlist->next ;
				vlist->next = nv ;
				ind[i] = numVerts ;
				pnode->index = numVerts ;
				numVerts ++ ;
			}

			if ( pnode->getSign( c1 ) == pnode->getSign( c2 ) )
			{
				sc[ i ] = 0 ;
			}
			else
			{
				sc[ i ] = 1 ;
			}
		}
		else
		{
			printf("Wrong!\n");
		}

	}

	if ( sc[ mini ] == 0 )
	{
//		printf("I am done!\n" ) ;
		return ;
	}

	/************************************************************************/
	/* Performing test                                                      */
	/************************************************************************/

	int fverts[4] ;
	int hasFverts[4] = { 0, 0, 0, 0 } ;
	int evert ;
	int needTess = 0 ;
	int location[4] ;
	int nvert[4] ={0,0,0,0};
	
	

	// First, face test
	int nbr[4][2] = { {0,1},{1,3},{2,3},{0,2} };
	int fdir[3][4] = {
		{2,1,2,1},
		{0,2,0,2},
		{1,0,1,0}};
	int dir3[3][4][2] = {
		{{1, -1},{2, 0},{1, 0},{2, -1}},
		{{2, -1},{0, 0},{2, 0},{0, -1}},
		{{0, -1},{1, 0},{0, 0},{1, -1}} };

	for ( i = 0 ; i < 4 ; i ++ )
	{
		int a = nbr[i][0];
		int b = nbr[i][1] ;

#ifndef TESS_UNIFORM
		if ( ht[a] != ht[b] )
#endif
		{
			// Different level, check if the dual edge passes through the face
			if ( hash->FindKey( (int) (node[a]), (int)(node[b]), fverts[i], location[i] ) )
			{
				// The vertex was found previously
				founds++ ;
				hasFverts[i] = 1 ;
				nvert[i] = 0 ;
				needTess = 1 ;
			}
			else
			{
				// Otherwise, we test it here
				int sht = ( ht[a] > ht[b] ? ht[b] : ht[a] ) ;
				int flen = ( 1 << sht ) ;
				int fst[3] ;

				fst[ fdir[dir][i] ] = st[ fdir[dir][i] ] ;
				fst[ dir3[dir][i][0] ] = st[ dir3[dir][i][0] ] + flen * dir3[dir][i][1] ;
				fst[ dir ] = st[ dir ] - ( st[ dir ] & (( 1 << sht ) - 1 ) ) ;

				if ( testFace( fst, flen, fdir[dir][i], mp[a], mp[b] ) == 0 )
				{
					// Dual edge does not pass face, let's make a new vertex
					VertexList* nv = new VertexList ;
					nv->vt[0] = 0 ;
					nv->vt[1] = 0 ;
					nv->vt[2] = 0 ;
					nv->next = vlist->next ;
					vlist->next = nv ;
					fverts[i] = numVerts ;
					location[i] = ((int) nv) ;
					nvert[i] = 1 ;
					numVerts ++ ;

					hash->InsertKey( (int)(node[a]), (int)(node[b]), fverts[i], location[i] ) ;


					hasFverts[ i ] = 1 ;
					needTess = 1 ;
					news ++ ;
				}
			}
		}
	}

	// Next, edge test
	int diag = 1 ;
	if ( needTess == 0 )
	{
		// Even if all dual edges pass through faces, the dual complex of an edge may not be convex
		//int st2[3] = { st[0], st[1], st[2] } ;
		//st2[ dir ] += len ;

		diag = testEdge( st, len, dir, node, mp ) ;
		if ( diag == 0 )
		{
			// When this happens, we need to create an extra vertex on the primal edge
			needTess = 1 ;
		}
	}

	float cent[3] ;
	if ( needTess )
	{
		edgeVerts ++ ;
		makeEdgeVertex( st, len, dir, node, mp, cent ) ;

		/* Just take centroid
		int num = 0 ;
		evert[0] = st[0] ;
		evert[1] = st[1] ;
		evert[2] = st[2] ;
		evert[dir] = 0 ;
		for ( i = 0 ; i < 4 ; i ++ )
		{
			evert[dir] += mp[i][dir] ;
			num ++ ;

			if ( hasFverts[ i ] )
			{
				evert[dir] += fverts[i][dir] ;
				num ++ ;
			}
		}
		evert[dir] /= num ;
		*/

/*
		if ( evert[dir] < st[dir] )
		{
			evert[dir] = st[dir] ;
		}
		else if ( evert[dir] > st[dir] + len )
		{
			evert[dir] = st[dir] + len ;
		}
*/		

		VertexList* nv = new VertexList ;
		nv->vt[0] = cent[0] ;
		nv->vt[1] = cent[1] ;
		nv->vt[2] = cent[2] ;
		nv->next = vlist->next ;
		vlist->next = nv ;
		evert = numVerts ;
		numVerts ++ ;

	}

	int flipped[3];
	if ( flip == 0 )
	{
		for (i = 0 ;i < 3 ; i ++)
		{
			flipped[i] = i ;
		}
	}
	else
	{
		for (i = 0 ;i < 3 ; i ++)
		{
			flipped[i] = 2 - i ;
		}
	}



	// Finally, let's output triangle
	if ( needTess == 0 )
	{
		// Normal splitting of quad
		if ( diag == 1 )
		{
			if ( node[0] != node[1] && node[1] != node[3] )
			{
				int tind1[]={0,1,3} ;
				
				numTris ++ ;
				IndexedTriangleList* t1 = new IndexedTriangleList ;
				t1->next = tlist->next;
				tlist->next = t1 ;
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					t1->vt[flipped[j]] = ind[ tind1[j] ] ;
				}
			}
			
			if ( node[3] != node[2] && node[2] != node[0] )
			{
				int tind2[]={3,2,0} ;
				
				numTris ++ ;
				IndexedTriangleList* t2 = new IndexedTriangleList ;
				t2->next = tlist->next;
				tlist->next = t2 ;
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					t2->vt[flipped[j]] = ind[ tind2[j] ] ;
				}
			}
		}
		else
		{
			if ( node[0] != node[1] && node[1] != node[2] )
			{
				int tind1[]={0,1,2} ;
				
				
				numTris ++ ;
				IndexedTriangleList* t1 = new IndexedTriangleList ;
				t1->next = tlist->next;
				tlist->next = t1 ;
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					t1->vt[flipped[j]] = ind[ tind1[j] ] ;
				}
			}
			
			if ( node[1] != node[3] && node[3] != node[2] )
			{
				int tind2[]={1,3,2} ;
				
				numTris ++ ;
				IndexedTriangleList* t2 = new IndexedTriangleList ;
				t2->next = tlist->next;
				tlist->next = t2 ;
				for ( int j = 0 ; j < 3 ; j ++ )
				{
					t2->vt[flipped[j]] = ind[ tind2[j] ] ;
				}
			}
		}
		
	}
	
	else
	{/*
		if ( flip == 1 )
		{
			int tempind[4]={ind[0], ind[2], ind[1], ind[3]};
			OctreeNode* tempnode[4] = {node[0], node[2], node[1], node[3]};
			int tempfverts[4]={fverts[3], fverts[2], fverts[1], fverts[0]};
			int temphasFverts[4]={hasFverts[3], hasFverts[2], hasFverts[1], hasFverts[0]};
			int templocation[4]={location[3], location[2], location[1], location[0]};
			int tempnvert[4]={nvert[3], nvert[2], nvert[1], nvert[0]};

			for ( i = 0 ; i < 4 ; i ++ )
			{
				ind[i] = tempind[i] ;
				node[i] = tempnode[i] ;
				fverts[i] = tempfverts[i] ;
				hasFverts[i] = temphasFverts[i] ;
				location[i] = templocation[i] ;
				nvert[i] = tempnvert[i] ;
			}
		}	
		*/
		int nnbr[4][2] = { {0,1},{1,3},{3,2},{2,0} };

		// Center-splitting
		for ( i = 0 ; i < 4 ; i ++ )
		{
			int a = nnbr[i][0];
			int b = nnbr[i][1] ;

			if ( hasFverts[ i ] )
			{
				// Further split each triangle into two
				numTris += 2 ;

				IndexedTriangleList* t = new IndexedTriangleList ;
				t->next = tlist->next;
				tlist->next = t ;
					t->vt[flipped[0]] = ind[ a ] ;
					t->vt[flipped[1]] = fverts[ i ] ;
					t->vt[flipped[2]] = evert ;

				t = new IndexedTriangleList ;
				t->next = tlist->next;
				tlist->next = t ;
					t->vt[flipped[0]] = evert ;
					t->vt[flipped[1]] = fverts[ i ] ;
					t->vt[flipped[2]] = ind[ b ] ;

				// Update geometric location of the face vertex
				VertexList* nv = ((VertexList *) location[i]) ;
				if ( nvert[i] )
				{
					nv->vt[0] = cent[0] ;
					nv->vt[1] = cent[1] ;
					nv->vt[2] = cent[2] ;
				}
				else
				{
					nv->vt[0] = ( nv->vt[0] + cent[0] ) / 2 ;
					nv->vt[1] = ( nv->vt[1] + cent[1] ) / 2 ;
					nv->vt[2] = ( nv->vt[2] + cent[2] ) / 2 ;
				}
			}
			else
			{
				// For one triangle with center vertex
				if ( node[a] != node[b] )
				{
					numTris ++ ;
					IndexedTriangleList* t = new IndexedTriangleList ;
					t->next = tlist->next;
					tlist->next = t ;
						t->vt[flipped[0]] = ind[ a ] ;
						t->vt[flipped[1]] = ind[ b ] ;
						t->vt[flipped[2]] = evert ;
				}
			}
		}
	}
//		printf("I am done!\n" ) ;
	
};
