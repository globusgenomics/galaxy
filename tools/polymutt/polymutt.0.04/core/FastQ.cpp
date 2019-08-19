#include "FastQ.h"

bool FastQ::OpenArchive(const char * filename)
   {
   CloseArchive();

   input = ifopen(filename, "rb");

   return input != NULL;
   }

void FastQ::CloseArchive()
   {
   if (input != NULL) ifclose(input);

   buffer.Clear();
   input = NULL;
   }

int FastQ::Load(IFILE & input, String & label, String & sequence, String & quality)
   {
   sequence.Clear();
   quality.Clear();

   if (buffer[0] == '>' || buffer[0] == '@')
      {
      label = buffer;
      buffer.Clear();
      }
   else
      label.ReadLine(input);

   // printf("LABEL: %s\n", (const char *) label);

   if (input == NULL || ifeof(input))
      return -1;

   if (label[0] != '>' && label[0] != '@')
      printf("WARNING: Invalid FAST[A/Q] sequence header (expecting line to begin with '>' or '@')\n");

   while (!ifeof(input))
      {
      buffer.ReadLine(input);
      // printf("BASES: %s\n", (const char *) buffer);

      if (buffer[0] == '+' || buffer[0] == '>' || buffer[0] == '@') break;

      int length = sequence.Length();
      for (int i = 0; i < buffer.Length(); i++)
         if (basekey[buffer[i]])
            {
            sequence += basekey[buffer[i]];
            length++;
            }

      sequence.SetLength(length);
      }

   if (label[0] == '@' && buffer[0] == '+')
      while (!ifeof(input) && quality.Length() < sequence.Length())
         {
         buffer.ReadLine(input);
         // printf("QUALS: %s\n", (const char *) buffer);

         // if (buffer[0] == '@') break;

         for (int i = 0; i < buffer.Length(); i++)
            if (buffer[i] > 32)
               quality += (char) (buffer[i] - char(32));
            else
               quality += (char) -1;

         quality.SetLength(buffer.Length());
         buffer.Clear();
         }

   return sequence.Length();
   }
