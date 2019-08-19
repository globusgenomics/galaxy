#ifndef __FAMILYLIKELIHOODSEQ_H__
#define __FAMILYLIKELIHOODSEQ_H__

#include "PedigreeGLF.h"
#include "NucFamGenotypeLikelihood.h"
#include "FamilyLikelihoodES.h"
#include "MutationModel.h"

//Def of FamGenotypeLikelihood class
class FamilyLikelihoodSeq : public NucFamGenotypeLikelihood
{
 public:
  FamilyLikelihoodES * fam;
  
 public:
  FamilyLikelihoodSeq();
  ~FamilyLikelihoodSeq();
  void   SetAlleleFreq(double);
  virtual double f(double freq);

  void InitFamilyLikelihoodES();
  void FillPenetrance();
  void FillPenetrance(FamilyLikelihoodES *famlk, PedigreeGLF *pedGLF);
  void FillZeroPenetrance(FamilyLikelihoodES *famlk, PedigreeGLF *pedGLF, int person, int genoIdx);
  double MonomorphismLogLikelihood_denovo(int, int);
  double PolymorphismLogLikelihood(int, int);
  double CalcAllFamLikelihood(double freq);
  double CalcAllFamLogLikelihood(double freq);
  double CalcSingleFamLikelihood(int i, double freq);
  double CalcSingleFamLikelihood_BA(int i, double freq);
  double CalcSingleFamLikelihood_denovo(int i, double freq);
  double CalcSingleFamLogLikelihood(int i, double freq);
  double CalcSingleFamLogLikelihood_BA(int i, double freq);
  double CalcSingleFamLogLikelihood_denovo(int i, double freq);
  void  CalcPostProb(double);
  void  CalcPostProb_SingleExtendedPed(int, double);
  void  CalcPostProb_SingleExtendedPed_BA(int, double);
  void  CalcPostProb_SingleExtendedPed_denovo(int, double);
};

#endif
