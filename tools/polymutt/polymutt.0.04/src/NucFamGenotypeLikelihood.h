#ifndef __NUCFAMGENOTYPELIKELIHOOD_H__
#define __NUCFAMGENOTYPELIKELIHOOD_H__

#include "CmdLinePar.h"
#include "PedigreeGLF.h"
#include "Parameters.h"
#include "StringArray.h"
#include "glfHandler.h"
#include "MathGold.h"
#include "MutationModel.h"
#include "Error.h"
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <vector>

//Def of NucFamGenotypeLikelihood class
class NucFamGenotypeLikelihood : public ScalarMinimizer
{
private:
  bool vcf;
public:
  int refBase;
  int allele1;
  int allele2;
  int geno11, geno12, geno22;
  int           currentFam;
  int           nFam;
  int           nPerson;
  int 	        nFounders;
  double        prior;
  double        priorFreq;
  double        theta;
  double        varPostProb;
  double        polyQual;
  int           totalDepth;
  double        avgDepth;
  int           numSampWithData;
  double        percSampWithData;
  double        avgMapQual;
  PedigreeGLF * pedGLF;
  Pedigree    * ped;
  double        parentPrior[9];
  double ***    kidsCondLikelihood;
  double **     transMissionProb;
  double **     parentMarginal;
  double **     parentConditional;
  double **     parentGLF;
  double ***    postProb;
  double **     dosage;
  int **        bestGenoIdx;
  String **     bestGenoLabel;
  bool          usePriorFreq;
  bool          useMleFreq;
  CmdLinePar *  par;
  std::vector<double> varllk;
  std::vector<double> varllk_noprior;
  std::vector<double> varfreq;
  AlleleMutationModel aM;
  GenotypeMutationModel gM;
  bool denovo_mono;
  double denovoLR;
  
public:
  NucFamGenotypeLikelihood();
  NucFamGenotypeLikelihood(PedigreeGLF *);
  ~NucFamGenotypeLikelihood();
  void   InitializeValues();
  void   InitializeTransmissionProb();
  void   InitializeParentMarginal();
  void   InitializePostProb();
  void   InitializeKidsCondLikelihood();
  void   DeleteTransmissionProb();
  void   DeletePostProb();
  void   DeleteParentMarginal();
  void   DeleteKidsCondLikelihood();
  void   SetCmdLinePar(CmdLinePar *p) { par = p; }
  int    GetFamCount();
  void   SetGLF(PedigreeGLF *);
  double f(double freq);
  void   SetAlleles(int, int);
  void   SetMleFreq(double freq); 
  void   SetParentPrior(double);
  void   SetParentPrior_denovo(double);
  void   SetParentPriorSingleTrio();
  void   SetParentPriorSingleTrio_denovo(double);
  void   SetTheta(double);
  void   SetPolyPrior();
  double GetPolyPrior();
  double GetPriorFreq();
  void   SetMleFreqFlag(bool);
  void   SetPriorFreqFlag(bool);
  void   SetPolyQual(double);
  void   SetDenovoMutationModel();
  void   CalcPolyQual(double postProb);
  void   CalcReadStats();
  double allFamLikelihood(double freq);
  double allFamLogLikelihood(double freq);
  double allFamLikelihood_denovo(double freq);
  double allFamLogLikelihood_denovo(double freq);
  double lkSingleFam(int i, double freq);
  double lkSingleFam_denovo(int i, double freq);
  double lkSinglePerson(int, int, double freq);
  double logLkSingleFam(int i, double freq);
  double logLkSingleFam_denovo(int i, double freq);
  double likelihoodKids(int, int, int, int, int); //conditional on parents
  double likelihoodONEKid(glfHandler *, int, int, int, int); //conditional on parents
  double likelihoodONEKid_denovo(glfHandler *, int, int, int, int); //conditional on parents allowing for de novo mutations
  double likelihoodKids_denovo(int, int, int, int, int); //conditional on parents allowing for de novo mutations
  double likelihoodKids_leaveone(int, int, int, int, int, int); //conditional on parents, leaving one kid and marginalizing others
  JointGenoLk KidJointGenoLikelihood(int, int);
  JointGenoLk_denovo KidJointGenoLikelihood_denovo(int, int);
  void   GetJointGenoLk_denovo(glfHandler *, int, int, int, int, double *);
  void   getLogGenoLikelihood(glfHandler * glf, int, int, unsigned char *, unsigned char *, unsigned char *);
  void   getGenoLikelihood(glfHandler * glf, int, int, double *, double *, double *);
  double OptimizeFrequency();
  double GetMaxLogLikelihood();
  double GetMaxLikelihood();
  double GetMinimizer();
  double PolymorphismLikelihood(int, int);
  double PolymorphismLogLikelihood(int, int);
  double PolymorphismLikelihoodSingleTrio(int, int);
  double MonomorphismLikelihood(int);
  double MonomorphismLogLikelihood(int);
  double MonomorphismLogLikelihood_denovo(int, int);
  double CalcDenovoMutLk(const double *, int, int);
  void   CalcParentMarginal(int, double);
  void   CalcParentMarginal_denovo(int, double);
  void   CalcParentMarginal_leaveone(int, double, int); // Marginalizing all kids except the ith
  void   CalcParentMarginalSingleTrio();
  void   CalcPostProb(double);
  void   CalcPostProb_SingleNucFam(int, double);
  void   CalcPostProb_SingleNucFam_denovo(int, double);
  void   CalcPostProb_SinglePerson(int, int, double);  
  JointGenoLk likelihoodKidGenotype(int, int, int, int, int, int);
  JointGenoLk_denovo likelihoodKidGenotype_denovo(int, int, int, int, int, int);
  void   GetBestGenotype(int, int, double);
  void   CalcDosage(int, int, double);
  int    GetBestGenoIdx(double p11, double p12, double p22);
  String GetBestGenoLabel(int best);
  String GetBestGenoLabel_denovo(int best);
  String GetBestGenoLabel_vcfv4(int best);
  double CalcDosage(int i,int j);
  void   toOriginalLikelihood(int loglk11, int loglk12, int loglk22, double *lk11, double *lk12, double *lk22);  double toOriginalLikelihood(int loglk); 
  int    CalcVarPosterior(int n);
  int    CalcMaxLogLkIdx(std::vector<double> &loglk, int n);
  double CalcSum(std::vector<double> &ratio);
  void   OutputVCF(FILE *);
  void   OutputVCF_denovo(FILE *);
};

#endif
