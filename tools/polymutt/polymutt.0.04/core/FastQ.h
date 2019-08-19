#ifndef __FASTQ__
#define __FASTQ__

#include "SequenceBasics.h"
#include "StringBasics.h"
#include "InputFile.h"

class FastQ : public SequenceBasics
   {
   public:
      IFILE    input;
      String   label;
      String   sequence;
      String   quality;

      FastQ()  { }
      ~FastQ() { CloseArchive(); }

      bool OpenArchive(const char * filename);
      void CloseArchive();

      int Load(IFILE & source, String & label, String & sequence, String & quality);
      int Load(IFILE & source)  { return Load(source, label, sequence, quality); }
      int Load()                { return Load(input); }

      int Length()     { return sequence.Length(); }

      bool isFastQ()   { return quality.Length() > 0; }
      bool isFastA()   { return !isFastQ(); }

      unsigned int GetWord(int start, int length)
         { return start + length <= Length() ? BitWord(&sequence[start], length) : 0; }

      unsigned int WordIsKnown(int start, int length)
         { return start + length <= Length() ? (FindAmbiguousBase(&sequence[start], length) < 0) : false; }

      String & PrintableSequence(String & output)
         { return SequenceBasics::PrintableSequence(output, &sequence[0], sequence.Length()); }

      String & PrintableQualities(String & output)
         { return SequenceBasics::PrintableQualities(output, &quality[0], sequence.Length()); }

   private:
      String buffer;
   };

#endif

