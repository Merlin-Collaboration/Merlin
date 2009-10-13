/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:55 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef deleters_h
#define deleters_h 1

#include <utility>

//	Function object for deleting containers of pointers.

template <class T>
class deleter 
{
  public:
      // Deletes the pointer p.
      void operator () (T* p);
};

//	Function object for deleting associative containers
//	whose value types are pointers.

template <class key, class val>
class map_deleter 
{
  public:
      void operator () (std::pair<key,val*>& arg);
};

template <class T>
inline void deleter<T>::operator () (T* p)
{
	if(p) delete p;
}

template <class key, class val>
inline void map_deleter<key,val>::operator () (std::pair<key,val*>& arg)
{
	if(arg.second) delete arg.second;
}

#endif
