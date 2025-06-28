/*
Description: A class cl2DArray of 2-dimensional matrices
*/
//==========================================================================

#ifndef CL2DARRAY_H
#define CL2DARRAY_H
#include <iostream>

//===========================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//============================================================== Declarations


//================================================================= cl2DArray
/* cl2DArray is a class of 2-dimensional arrays. It generalizes (but not
 * inherits) clArray in most aspects. cl2DArray is mainly introduced to
 * realize 2-dimensional arrays in the least CPU-time consuming way.
 */

template <class T>
class cl2DArray
{private:
  T *arr;                       //pointer to the array with the values
  long unsigned int lenarr;     //cumulative size of the array with the values


public:
  int  First1, First2, Last1, Last2;          //1st and last indices
  unsigned int N1, N2;          //Sizes of the array
cl2DArray(): N1(1), N2(1), First1(0), First2(0), Last1(0), Last2(0)
    {};

  cl2DArray(const unsigned int N1in, const unsigned int N2in);
  int SetFstI(const int First1in=0, const int First2in=0);
  int Reset(const unsigned int N1in=1, const unsigned int N2in=1);
  T& operator()(const int i1, const int i2);

  ~cl2DArray(){
      delete [] arr;
    } //~cl2DArray -> 
  
};//cl2DArray ->



//================================================= cl2DArray
template <class T>
cl2DArray<T>::cl2DArray(const unsigned int N1in, const unsigned int N2in)
{
  N1=N1in; N2=N2in; First1=First2=0; Last1=N1-1; Last2=N2-1; 
  lenarr=N1*N2;
  //Allocate the memory for the array
  // with the values
  arr = new T[lenarr];
  
}//cl2DArray<T>::cl2DArray ->

//--------------------------------------------------------------------------
template <class T>
inline int cl2DArray<T>::SetFstI(const int First1in, const int First2in)
// Sets the values of first indices
{
  First1 = First1in; First2 = First2in; Last1 = First1 +N1 -1; Last2 = First2 +N2 -1;
  return 0;
}//cl2Darray<T>::SetFstI ->

//--------------------------------------------------------------------------
template <class T>
inline int cl2DArray<T>::Reset(const unsigned int N1in,
                               const unsigned int N2in)
{
  delete [] arr;

  N1 = N1in; N2 = N2in; First1 = First2 =0; Last1 = N1-1; Last2 = N2-1;
  lenarr = N1*N2;
  arr = new T[lenarr];

  return 0;
}//cl2DArray<T>::Reset->

//---------------------------------------------------------------------------
template <class T>
T& cl2DArray<T>::operator()(const int k1,  const int k2)
// Returns the value of the array with given indices
{ long int offset;
  //Calculate the offset in arr
  offset = N2*(k1 - First1) + k2-First2;
  if (offset < 0 || offset > (long int)lenarr-1)                        //Prevent overflow
    {
      std::cout << "cl2DArray: Wrong indices "<<k1<<", "<<k2<<" in cl2DArray"<<std::endl;
      std::cout << "Allowed ranges of indices are "<<First1<<" < k1 < "<<Last1<<" "<<First2<<" < k2 < "<<Last2<<std::endl; 
    }

  return *(arr+offset);
}//cl2DArray<T>::operator() ->

#endif //CL2DARRAY_H
