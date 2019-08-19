#ifndef __INPUTFILE_H__
#define __INPUTFILE_H__

#ifdef  __gnu_linux__
#ifndef __ZLIB_AVAILABLE__
#define __ZLIB_AVAILABLE__
#endif
#endif

#ifdef  __ZLIB_AVAILABLE__

#include <zlib.h>
#include <stdio.h>

class IFILE
   {
   public:
      bool gzMode;
      union
         {
         gzFile gzHandle;
         FILE * handle;
         };

   IFILE()
      {
      gzMode = false;
      handle = NULL;
      }

   IFILE(const char * filename, const char * mode);

   operator void * ()
      { return gzMode ? (void *) gzHandle : (void *) handle; }

   IFILE operator = (const IFILE & rhs)
      {
      if ((gzMode = rhs.gzMode) == true)
         gzHandle = rhs.gzHandle;
      else
         handle = rhs.handle;

      return *this;
      }

   IFILE operator = (FILE * rhs)
      {
      gzMode = false;
      handle = rhs;
      return *this;
      }

   IFILE operator = (gzFile & rhs)
      {
      gzMode = true;
      gzHandle = rhs;
      return *this;
      }

   bool operator == (void * rhs)
      {
      if (rhs != NULL)
         return false;
      return gzMode ? gzHandle == rhs : handle == rhs;
      }
   };

inline IFILE ifopen(const char * filename, const char * mode)
   { IFILE file(filename, mode); return file; }

inline int ifclose(IFILE & file)
   {
   int result = file.gzMode ? gzclose(file.gzHandle) : fclose(file.handle);
   file.gzHandle = NULL;
   return result;
   }

inline int ifgetc(IFILE & file)
   { return file.gzMode ? gzgetc(file.gzHandle) : fgetc(file.handle); }

inline void ifrewind(IFILE & file)
   { if (file.gzMode) gzrewind(file.gzHandle); else rewind(file.handle); }

inline int ifeof(IFILE & file)
   { return file.gzMode ? gzeof(file.gzHandle) : feof(file.handle); }

inline unsigned int ifread(IFILE & file, void * buffer, unsigned int size)
   { return file.gzMode ? gzread(file.gzHandle, buffer, size) :
                          fread(buffer, 1, size, file.handle); }

inline unsigned int ifwrite(IFILE & file, void * buffer, unsigned int size)
   { return file.gzMode ? gzwrite(file.gzHandle, buffer, size) :
                          fwrite(buffer, 1, size, file.handle); }

#else

#include <stdio.h>

class IFILE
   {
   public:
      FILE * handle;

      IFILE()
         { handle = NULL; }
      IFILE(const char * filename, const char * mode)
         { handle = fopen(filename, mode); }
      ~IFILE()
         { }

      operator FILE *()
         { return handle; }

      IFILE & operator = (FILE * rhs)
         { handle = rhs; return *this; }

      IFILE & operator = (const IFILE & rhs)
         { handle = rhs.handle; return * this; }

      bool operator == (void * rhs)
         {
         if (rhs != NULL)
            return false;
         return handle == rhs;
         }
   };

inline IFILE ifopen(const char * filename, const char * mode)
   { IFILE file(filename, mode); return file; }

inline int ifclose(IFILE & file)
   {
   int result = fclose(file.handle);
   file.handle = NULL;
   return result;
   }

inline int ifgetc(IFILE & file)
   { return fgetc(file.handle); }

inline void ifrewind(IFILE & file)
   { rewind(file.handle); }

inline int ifeof(IFILE & file)
   { return feof(file.handle); }

inline unsigned int ifread(IFILE & file, void * buffer, unsigned int size)
   { return fread(buffer, 1, size, file.handle); }

inline unsigned int ifwrite(IFILE & file, void * buffer, unsigned int size)
   { return fwrite(buffer, 1, size, file.handle); }

#endif

int ifprintf(IFILE output, char * format, ...);

#endif

