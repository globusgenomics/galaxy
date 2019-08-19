#include "CustomSettings.h"

#define   __CS_INTEGER__     0
#define   __CS_DOUBLE__      1
#define   __CS_STRING__      2

class CustomSetting
   {
   public:
      CustomSetting(int t, void * v)
         {type = t; value = v; }

      int    type;
      void * value;
   };

CustomSettings::CustomSettings()
   {

   }

CustomSettings::~CustomSettings()
   {
   for (int i = 0; i < settings.Capacity(); i++)
      if (settings.SlotInUse(i))
         delete (CustomSetting *) settings.Object(i);
   }

void CustomSettings::ShowHeader()
   {
   if (settings.Entries() == 0)
      {
      printf("The default settings are ...\n");
      }
   }

void CustomSettings::AddSetting(const char * name, double & value)
   {
   String string(name);

   if (settings.Find(string) >= 0)
      delete (CustomSetting *) settings.Object(string);

   ShowHeader();

   settings.SetObject(string, new CustomSetting(__CS_DOUBLE__, &value));

   printf("  %s is set to %.3g\n", name, value);
   }

void CustomSettings::AddSetting(const char * name, int & value)
   {
   String string(name);

   if (settings.Find(string) >= 0)
      delete (CustomSetting *) settings.Object(string);

   ShowHeader();

   settings.SetObject(string, new CustomSetting(__CS_INTEGER__, &value));

   printf("  %s is set to %d\n", name, value);
   }

void CustomSettings::AddSetting(const char * name, String & value)
   {
   String string(name);

   if (settings.Find(string) >= 0)
      delete (CustomSetting *) settings.Object(string);

   ShowHeader();

   settings.SetObject(string, new CustomSetting(__CS_STRING__, &value));

   printf("  %s is set to %s\n", name, (const char *)value);
   }

void CustomSettings::LoadSettings(const char * filename)
   {
   FILE * file = fopen(filename, "rt");

   if (file == NULL)
      {
      printf("Settings file not available, using default settings\n\n");
      return;
      }

   printf("Loading custom settings from file %s ...\n", filename);

   LoadSettings(file);

   fclose(file);
   }

void CustomSettings::LoadSettings(FILE * file)
   {
   String input;
   String variable;
   String value;

   while (!feof(file))
      {
      input.ReadLine(file);
      input.Trim();

      if (input.Length() == 0 || input[0] == '#' || input[0] == '/' && input[1] == '/')
         continue;

      int equal_sign = input.FindChar('=');

      if (equal_sign < 0)
         {
         printf("  INCORRECT INPUT: %s\n", (const char *) input);
         continue;
         }

      int setting = settings.Find(input.Left(equal_sign).Trim());

      if (setting < 0)
         {
         printf("  UNKNOWN SETTING: %s\n", (const char *) input);
         continue;
         }

      input = input.SubStr(equal_sign + 1).Trim();

      CustomSetting * ptr = (CustomSetting *) settings.Object(setting);

      switch (ptr->type)
         {
         case __CS_INTEGER__ :
            * (int *) ptr->value = input.AsInteger();
            printf("  Setting %s changed to %d\n",
                   (const char *) settings[setting], * (int *) ptr->value);
            break;
         case __CS_DOUBLE__ :
            * (double *) ptr->value = input.AsDouble();
            printf("  Setting %s changed to %.3g\n",
                   (const char *) settings[setting], * (double *) ptr->value);
            break;
         case __CS_STRING__ :
	    * (String *) ptr->value = input;
	    printf("  Setting %s changed to %s\n", 
                    (const char *) settings[setting], 
                    (const char *)(* (String *) ptr->value));
         default:
         printf("    INTERNAL ERROR: %s\n", (const char *) input);
         }
      }
   printf("\n");
   }

   
