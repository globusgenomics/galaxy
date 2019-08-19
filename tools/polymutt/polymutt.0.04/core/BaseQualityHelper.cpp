#include "BaseQualityHelper.h"

#include <math.h>

baseQualityConvertor bQualityConvertor;

baseQualityConvertor::baseQualityConvertor()
   {
   // Create a quick lookup table to speed up conversion of
   // base quality values stored as log10 (error rates) into
   // fractional error rates
   for (int i = 0; i < 255; i++)
      doubleLookup[i] = pow(0.1, i * 0.1);
   doubleLookup[255] = 0.0;
   }

double baseQualityConvertor::toDouble(unsigned char bq)
   {
   return doubleLookup[bq];
   }

