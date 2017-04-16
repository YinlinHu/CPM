#pragma once
#include "stdio.h"
#include <vector>

template <class T>
void _Release1DBuffer(T* pBuffer)
{
	if(pBuffer!=NULL)
		delete []pBuffer;
	pBuffer=NULL;
}

template <class T>
void _Rlease2DBuffer(T** pBuffer,size_t nElements)
{
	for(size_t i=0;i<nElements;i++)
		delete [](pBuffer[i]);
	delete []pBuffer;
	pBuffer=NULL;
}

#ifdef _MATLAB
#include "mex.h"
#endif


#ifndef WIN32

#define strcmpi strcasecmp

template <class T1,class T2>
T1 __min(T1 a, T2 b)
{
  return (a>b)?b:a;
}

template <class T1,class T2>
T1 __max(T1 a, T2 b)
{
  return (a<b)?b:a;
}

#endif
