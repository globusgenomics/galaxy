#ifndef __BASEQUALITY_H__
#define __BASEQUALITY_H__

class baseQualityConvertor
   {
   public:
      baseQualityConvertor();

      double toDouble(unsigned char baseQuality);

   private:
      double doubleLookup[256];
   };

extern baseQualityConvertor bQualityConvertor;


#endif


