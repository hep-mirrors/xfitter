#ifndef FANTOMASPRINT_H
#define FANTOMASPRINT_H

extern double fantomuR;
extern double fantomuF;

inline void trackLowestValue(double *newValue, double &lowestValue){
  if (*newValue < lowestValue) 
    lowestValue = *newValue;
}
#endif
