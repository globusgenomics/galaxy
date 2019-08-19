#ifndef __FamilyLikelihoodES_H__
#define __FamilyLikelihoodES_H__

#include "MutationModel.h"
#include "glfHandler.h"
#include "Pedigree.h"
#include "MathMatrix.h"
#include <vector>
#include <map>

using namespace std;


class ES_Peeling
{
public:
 Pedigree *ped;
 Family *family;
 int famIdx;
 int famSize;
 IntArray *parents;
 IntArray *offspring;
 IntArray *spouses;
 IntArray leaf;
 std::vector<std::pair<int, int> > roof;
 IntArray peripheral;
 std::vector<std::pair<int, int> > from;
 std::vector<std::pair<int, int> > to;
 std::vector<int> peelingType; //1: offspring to parents, 2: one couple to another, 3: parents to child
  
 public:
 ES_Peeling();
 ~ES_Peeling();
 
 void SetPedigree(Pedigree *);
 void SetFamily(Family *);
 void SetFamily(int);
 void SetupConnections();
 void BuildInitialPeelable();
 void BuildPeelingOrder();
 void BuildInitialPeelable2();
 void BuildPeelingOrder2();
 void PrintPeelingOrder();
 void PrintPeelable();
 int find_element(std::vector<pair<int, int> >&, std::pair<int, int>&);
 bool remove_element(std::vector<pair<int, int> >&, int);
 bool UpdateRoof(std::vector<pair<int, int> > &roof, int index);
 bool isRoof(int);
 bool isLeaf(int);
 bool isPeripheral(int);
 bool isFinal(int);
 String GetPID(int);      
};


class FamilyLikelihoodES
{
 public:
  Pedigree *ped;
  Family *family;
  int famIdx;
  int famSize;
  int nFounders;
  int allele1, allele2;
  double frequency;
  std::vector<int> genoIdx;

  //de novo mutation models
  AlleleMutationModel *aM;
  GenotypeMutationModel *gM;
      
 ES_Peeling es;
 
 vector<double> priors;
 Matrix penetrances;
 Matrix states; //indices of genotypes
 vector<double> ** transmission;
 vector<double> ** transmission_denovo;
 vector<double> ** transmission_BA;
 Matrix partials;
 std::map<pair<int, int>, vector<vector<double> > > marriage_partials;
             
 public:
 FamilyLikelihoodES();
 ~FamilyLikelihoodES();
 
 void InitValues();
 void SetPedigree(Pedigree *);
 void SetFamily(Family *);
 void SetFamilyIndex(int);
 void PreparePeeling();
 bool isPhenotyped(int);
 void SetAlleles(int a1, int a2);
 void SetFounderPriors(double);
 void SetFounderPriors_BA(double);
 void SetTransmissionMatrix();
 void SetTransmissionMatrix_denovo();
 void SetTransmissionMatrix_BA();
 void SetMarriagePartials(std::pair<int, int>&);
 void SetMarriagePartials_BA(std::pair<int, int>&);
 void SetGenotypeMutationModel(GenotypeMutationModel *);
 void PrintTransmissionMatrix();
 void PrintTransmissionMatrix_BA();
 void FreeTransmissionMatrix(); 
 void FreeTransmissionMatrix_denovo(); 
 void FreeTransmissionMatrix_BA(); 
 void InitializeStates(double freq);
 void InitializePartials();
 void InitializePartials_BA();
 void FillPenetrance(Family *);
 void FillPenetrance();
 void FillPenetrance(int);
 void PrintPenetrance();
 void PrintPenetrance_BA();
 double CalculateLikelihood();
 double CalculateLikelihood_BA();
 double CalculateLikelihood_denovo();
 void CalculateLikelihood(Family *);
 void CalculateLikelihood(int);
 void peelParents2Offspring(int);
 void peelOffspring2Parents(int);
 void peelSpouse2Spouse(int);
 void peelParents2Offspring_BA(int);
 void peelOffspring2Parents_BA(int);
 void peelSpouse2Spouse_BA(int);
 void peelParents2Offspring_denovo(int);
 void peelOffspring2Parents_denovo(int);
 void peelSpouse2Spouse_denovo(int);
 String GetPID(int);
};

#endif
