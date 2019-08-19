#include "SequenceBasics.h"

char SequenceBasics::basekey[256] = {
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,   0,  -2,  -2,  -2,  -2,  -2,  -2,  /* 0 -> 16 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 16 -> 32 */
       0,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 32 -> 48 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 48 -> 64 */
      -2,   1,  -2,   2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -1,  -2,  /* 64 -> 80 */
      -2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 80 -> 96 */
      -2,   1,  -2,   2,  -2,  -2,  -2,   3,  -2,  -2,  -2,  -2,  -2,  -2,  -1,  -2,  /* 96 -> 112 */
      -2,  -2,  -2,  -2,   4,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 112 -> 128 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 128 -> 144 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 144 -> 160 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 160 -> 176 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 176 -> 192 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 192 -> 208 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 208 -> 224 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  /* 224 -> 240 */
      -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2   /* 240 -> 256 */
     };

unsigned int SequenceBasics::BitWord(const char * sequence, int length)
   {
   unsigned int hash = 0;

   for (int i = 0; i < length; i++, sequence++)
      hash = hash * 4 + *sequence - 1;

   return hash;
   }

String & SequenceBasics::TranslateBitWord(String & string, int bitWord, int length)
   {
   char base[] = {'A', 'C', 'G', 'T'};

   string.Clear();

   for (int i = 0; i < length; i++, bitWord /= 4)
      string += (char) base[bitWord & 3];
   string.Reverse();

   return string;
   }

String & SequenceBasics::PrintableSequence(String & string, const char * sequence, int length)
   {
   char base[] = {'N', 'A', 'C', 'G', 'T'};

   string.Clear();

   for (int i = 0; i < length; i++, sequence++)
      string += (char) base[*sequence > 0 ? *sequence : 0];

   return string;
   }

String & SequenceBasics::PrintableQualities(String & string, const char * quality, int length)
   {
   string.Clear();

   for (int i = 0; i < length; i++, quality++)
      string += (char) (*quality + 32);

   return string;
   }

int SequenceBasics::FindAmbiguousBase(const char * sequence, int length)
   {
   for (int i = length - 1; i >= 0; i--)
      if (sequence[i] < 0)
         return i;

   return -1;
   }

