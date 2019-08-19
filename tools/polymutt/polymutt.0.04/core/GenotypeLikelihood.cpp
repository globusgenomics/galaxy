#include "BaseQualityHelper.h"
#include "StringArray.h"
#include "Parameters.h"
#include "glfHandler.h"
#include "MathGold.h"
#include "Error.h"

#include <math.h>
#include <time.h>

class GenotypeLikelihood : public ScalarMinimizer
   {
   public:
      int n;
      glfHandler * glf;

      void SetAlleles(int al1, int al2)
         {
         allele1 = al1;
         allele2 = al2;

         geno11 = glfHandler::GenotypeIndex(allele1, allele1);
         geno12 = glfHandler::GenotypeIndex(allele1, allele2);
         geno22 = glfHandler::GenotypeIndex(allele2, allele2);
         }

      virtual double f(double freq) { return -Evaluate(freq); }

      double Evaluate(double freq)
         {
         double prior11 = freq * freq;
         double prior12 = freq * (1.0 - freq) * 2.0;
         double prior22 = (1.0 - freq) * (1.0 - freq);

         double likelihood = 1.0;

         for (int i = 0; i < n; i++)
            likelihood *= prior11 * glf[i].likelihoods[geno11] +
                          prior12 * glf[i].likelihoods[geno12] +
                          prior22 * glf[i].likelihoods[geno22];

         return likelihood;
         }

      void GetPriors(double * priors, double freq)
         {
         for (int i = 0; i < 10; i++)
            priors[i] = 0.0;

         priors[geno11] = freq * freq;
         priors[geno12] = freq * (1.0  - freq) * 2.0;
         priors[geno22] = (1.0 - freq) * (1.0 - freq);
         }

      double OptimizeFrequency()
         {
         a = 0.00001; fa = f(a);
         b = 0.99999; fb = f(b);
         c = 0.5; fc = f(c);

         Brent(0.0001);

         return min;
         }

   protected:
      int allele1, allele2;
      int geno11, geno12, geno22;
   };

FILE * baseCalls = NULL;

void DumpDetails(glfHandler * glf, int n)
   {
   char alleles[] = { 0, 'a', 'c', 'g', 't' };

   //printf("Dump for section %s, position %d [%c]\n",
   //       (const char *) glf[0].label, glf[0].currentEntry, alleles[glf[0].data.refBase]);

   printf("Depth");
   for (int i = 0; i < n; i++)
      printf("\t%d", glf[i].data.depth);
   printf("\n");

   printf("MapQ");
   for (int i = 0; i < n; i++)
      printf("\t%d", glf[i].data.mapQuality);
   printf("\n");

   for (int i = 1, index = 0; i <= 4; i++)
      for (int j = i; j <= 4; j++, index++)
         {
         printf("%c/%c", alleles[i], alleles[j]);
         for (int k = 0; k < n; k++)
            printf("\t%d", glf[k].data.lk[index]);
         printf("\n");
         }
   }

int GetBestGenotype(double likelihoods[], double priors[])
   {
   int best = 0;

   for (int i = 1; i < 10; i++)
      if (likelihoods[i] * priors[i] > likelihoods[best] * priors[best])
         best = i;

   return best;
   }

const char * GetBestGenotypeLabel(double likelihoods[], double priors[])
   {
   const char * genotypeLabel[10] = {"A/A", "A/C", "A/G", "A/T", "C/C", "C/G", "C/T", "G/G", "G/T", "T/T"};

   return genotypeLabel[GetBestGenotype(likelihoods, priors)];
   }

int GetBestRatio(double likelihoods[], double priors[])
   {
   double sum = 0.0;
   int best = 0;

   for (int i = 1; i < 10; i++)
      if (likelihoods[i] * priors[i] > likelihoods[best] * priors[best])
         best = i;

   for (int i = 0; i < 10; i++)
      sum += likelihoods[i] * priors[i];

   double error = 1.0 - likelihoods[best] * priors[best]/sum;

   if (error < 0.0000000001)
      return 100;

   return int (-log10(error) * 10 + 0.5);
   }

void ReportGenotypes(glfHandler * glf, int n, int al1, int al2)
   {
   if (baseCalls == NULL)
      return;

   double priors[10];
   GenotypeLikelihood lk;

   // Find best frequency
   lk.glf = glf;
   lk.n = n;
   lk.SetAlleles(al1, al2);
   lk.OptimizeFrequency();
   lk.GetPriors(priors, lk.min);

   fprintf(baseCalls, "\t%.3f", lk.min);

   for (int i = 0; i < n; i++)
      {
      // Report on the best genotype for the current SNP model
      int quality = GetBestRatio(glf[i].likelihoods, priors);

      fprintf(baseCalls, "\t%s\t%d", GetBestGenotypeLabel(glf[i].likelihoods, priors), quality);
      }

   fprintf(baseCalls, "\n");
   }

void ReportLocation(glfHandler * glf, int totalCoverage, double averageMapQuality)
   {
   char alleles[] = { 0, 'a', 'c', 'g', 't' };

  // if (baseCalls != NULL)
      //fprintf(baseCalls, "%s\t%d\t%c\t%d\t%.2f\t",
      //       (const char *) glf[0].label, glf[0].currentEntry,
      //       alleles[glf[0].data.refBase],
      //       totalCoverage, averageMapQuality);
}

void ReportSNP(glfHandler * glf, int n, int allele1, int allele2, double posterior)
   {
   char alleles[] = { 0, 'a', 'c', 'g', 't' };

   int quality = posterior > 0.9999999999 ? 100 : int (-log10(1 - posterior) * 10 + 0.5);

   if (baseCalls != NULL)
      fprintf(baseCalls, "%c/%c\t%d", alleles[allele1], alleles[allele2], quality);

   ReportGenotypes(glf, n, allele1, allele2);
   }

double PolymorphismLikelihood
         (glfHandler * glf, int n, int refAllele, int mutAllele)
   {
   GenotypeLikelihood lk;

   lk.glf = glf;
   lk.n = n;
   lk.SetAlleles(refAllele, mutAllele);
   lk.OptimizeFrequency();

   return lk.fmin;
   }

double SinkLikelihood
       (glfHandler * glf, int n)
   {
   double lk = 0.0;

   for  (int r = 1; r <= 4; r++)
      for (int m = r + 1; m <= 4; m++)
          {
          int geno = glfHandler::GenotypeIndex(r, m);

          double partial = 1.0;
          for (int i = 0; i < n; i++)
             partial *= glf[i].likelihoods[geno];

          lk += partial;
          }
 
   return lk / 6.;
   }

