#ifndef __CUSTOMSETTINGS_H__
#define __CUSTOMSETTINGS_H__

#include "StringHash.h"

class CustomSettings
   {
   public:
      CustomSettings();
      ~CustomSettings();

      void ShowHeader();
      void AddSetting(const char * label, double & value);
      void AddSetting(const char * label, int & value);
      void AddSetting(const char * label, String & value);

      void LoadSettings(const char * filename);
      void LoadSettings(FILE * file);

      void ListSettings();

   private:
      StringHash settings;
   };


#endif

