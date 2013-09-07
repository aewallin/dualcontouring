/*

  A hash table hashed by two/four integers.

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


#ifndef HASHMAP_H
#define HASHMAP_H

#include <stdio.h>
#include <stdlib.h>

#define HASH_LENGTH  20
#define MAX_HASH (1<<HASH_LENGTH)
#define HASH_BIT_2 10
#define HASH_BIT_4 5

struct ElementList
{
	int n1;
	int n2;
	int index ;
	int location ;
	ElementList* next;
};

struct ElementList2
{
	int n1;
	int n2;
	float v[3] ;
	ElementList2* next;
};

struct ElementList4
{
	int n1;
	int n2;
	int n3;
	int n4;
	int index ;
	ElementList4* next;
};

class HashMap
{
/// Hash table
	ElementList *table[MAX_HASH];

/// Create hash key
	int CreateKey( int k1, int k2 )
	{
		int ind = (	(( k1 & ( (1<<HASH_BIT_2) - 1 )) << ( HASH_BIT_2 )) |
					( k2 & ( (1<<HASH_BIT_2) - 1)) ) 
					& ((1<<HASH_LENGTH) - 1);

		return ind;
	}

public:

	/// Constructor
	HashMap ( )
	{
		for ( int i = 0; i < MAX_HASH; i ++ )
		{
			table[i] = NULL;
		}
	};

	/// Lookup Method
	int FindKey( int k1, int k2, int& index, int& location )
	{
		/// Create hash key
		int ind = CreateKey ( k1, k2 );

		/// Find it in the table
		ElementList *p = table[ind];
		while (p)
		{
			if ((p->n1 == k1) && (p->n2 == k2))
			{
				index = p->index ;
				location = p->location ;
				return 1 ;
			}
			p = p->next;
		}

		return 0 ;
	};


	/// Insertion method
	void InsertKey ( int k1, int k2, int index, int location )
	{
		/// Create hash key
		int ind = CreateKey ( k1, k2 );

		/// Create hash entry
		ElementList *node = new ElementList;
		
		node->index = index ;
		node->location = location ;
		node->n1 = k1 ;
		node->n2 = k2 ;
		node->next = table[ ind ] ;
		table[ ind ] = node ;
	};

	// Destruction method
	~HashMap()
	{
		ElementList *p, *pp;

		for ( int i = 0; i < MAX_HASH; i ++)
		{
			p = table[i];

			while (p)
			{
				pp = p->next;
				delete p;
				p = pp;
			}

		}

	};

};

class HashMap2
{
/// Hash table
	ElementList2 *table[MAX_HASH];

/// Create hash key
	int CreateKey( int k1, int k2 )
	{
		int ind = (	(( k1 & ( (1<<HASH_BIT_2) - 1 )) << ( HASH_BIT_2 )) |
					( k2 & ( (1<<HASH_BIT_2) - 1)) ) 
					& ((1<<HASH_LENGTH) - 1);

		return ind;
	}

public:

	/// Constructor
	HashMap2 ( )
	{
		for ( int i = 0; i < MAX_HASH; i ++ )
		{
			table[i] = NULL;
		}
	};

	/// Lookup Method
	int FindKey( int k1, int k2, float coord[3] )
	{
		/// Create hash key
		int ind = CreateKey ( k1, k2 );

		/// Find it in the table
		ElementList2 *p = table[ind];
		while (p)
		{
			if ((p->n1 == k1) && (p->n2 == k2))
			{
				coord[0] = p->v[0] ;
				coord[1] = p->v[1] ;
				coord[2] = p->v[2] ;
				return 1 ;
			}
			p = p->next;
		}

		return 0 ;
	};


	/// Insertion method
	void InsertKey ( int k1, int k2, float coord[3] )
	{
		/// Create hash key
		int ind = CreateKey ( k1, k2 );

		/// Create hash entry
		ElementList2 *node = new ElementList2;
		
		node->v[0] = coord[0] ;
		node->v[1] = coord[1] ;
		node->v[2] = coord[2] ;
		node->n1 = k1 ;
		node->n2 = k2 ;
		node->next = table[ ind ] ;
		table[ ind ] = node ;
	};

	// Destruction method
	~HashMap2()
	{
		ElementList2 *p, *pp;

		for ( int i = 0; i < MAX_HASH; i ++)
		{
			p = table[i];

			while (p)
			{
				pp = p->next;
				delete p;
				p = pp;
			}

		}

	};

};

class HashMap4
{
/// Hash table
	ElementList4 *table[MAX_HASH];

/// Create hash key
	int CreateKey( int k1, int k2, int k3, int k4 )
	{
		int ind = (	(( k1 & ( (1<<HASH_BIT_4) - 1 )) << ( 3 * HASH_BIT_4 )) |
					(( k2 & ( (1<<HASH_BIT_4) - 1 )) << ( 2 * HASH_BIT_4 )) |
					(( k3 & ( (1<<HASH_BIT_4) - 1 )) << ( HASH_BIT_4 )) |
					( k4 & ( (1<<HASH_BIT_4) - 1)) ) 
					& ((1<<HASH_LENGTH) - 1);

		return ind;
	}

public:

	/// Constructor
	HashMap4 ( )
	{
		for ( int i = 0; i < MAX_HASH; i ++ )
		{
			table[i] = NULL;
		}
	};

	/// Lookup Method
	int FindKey( int k1, int k2, int k3, int k4 )
	{
		/// Create hash key
		int ind = CreateKey ( k1, k2, k3, k4 );

		/// Find it in the table
		ElementList4 *p = table[ind];
		while (p)
		{
			if ((p->n1 == k1) && (p->n2 == k2) && (p->n3 == k3) && (p->n4 == k4))
			{
				return p->index ;
			}
			p = p->next;
		}

		return -1 ;
	};


	/// Insertion method
	void InsertKey ( int k1, int k2, int k3, int k4, int index )
	{
		/// Create hash key
		int ind = CreateKey ( k1, k2, k3, k4 );

		/// Create hash entry
		ElementList4 *node = new ElementList4;
		
		node->index = index ;
		node->n1 = k1 ;
		node->n2 = k2 ;
		node->n3 = k3 ;
		node->n4 = k4 ;
		node->next = table[ ind ] ;
		table[ ind ] = node ;
	};

	// Destruction method
	~HashMap4()
	{
		ElementList4 *p, *pp;

		for ( int i = 0; i < MAX_HASH; i ++)
		{
			p = table[i];

			while (p)
			{
				pp = p->next;
				delete p;
				p = pp;
			}

		}

	};

};

#endif