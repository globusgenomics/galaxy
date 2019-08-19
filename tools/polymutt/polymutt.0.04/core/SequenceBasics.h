#ifndef __SEQUENCEBASICS_H__
#define __SEQUENCEBASICS_H__

#include "StringBasics.h"

class SequenceBasics
   {
   public:
      static char basekey[256];

      static unsigned int BitWord(const char * sequence, int length);
      static String & TranslateBitWord(String & string, int bitWord, int length);
      
      static String & PrintableSequence(String & string, const char * sequence, int length);
      static String & PrintableQualities(String & string, const char * quality, int length);

      static int FindAmbiguousBase(const char * sequence, int length);

   };

#endif

