#ifndef __PEDIGREEGLF_H__
#define __PEDIGREEGLF_H__

#include "glfHandler.h"
#include "Pedigree.h"
#include "BaseQualityHelper.h"
#include "StringMap.h"
#include "Error.h"
#include <vector>

class Poly
{
public:
  // Create convenient aliases for each base
  static int ts(int refBase) {
    //unsigned char ts = (((refBase - 1) ^ 2) + 1);
    int ts = 0; //= (((refBase - 1) ^ 2) + 1);

    switch(refBase) {
    case 1:  ts = 3; break;
    case 2:  ts = 4; break;
    case 3:  ts = 1; break;
    case 4:  ts = 2; break;
    }
    //    printf("in ts() ref==%d ts=%d\n", refBase, ts);
    return(ts);
  }
  static int tvs1(int refBase) {

    //unsigned tvs1 = (((refBase - 1) ^ 3) + 1);
    int tvs1 = 0; // = (((refBase - 1) ^ 3) + 1);
    switch(refBase) {
    case 1:  tvs1 = 2; break;
    case 2:  tvs1 = 1; break;
    case 3:  tvs1 = 2; break;
    case 4:  tvs1 = 1; break;
    }
    //    printf("in tvs1() ref==%d tvs1=%d\n", refBase, tvs1);
    return(tvs1);
  }
  static int tvs2(int refBase) {

    //unsigned char tvs2 = (((refBase - 1) ^ 1) + 1);
    int tvs2 = 0 ; //= (((refBase - 1) ^ 1) + 1);
    switch(refBase) {
    case 1:  tvs2 = 4; break;
    case 2:  tvs2 = 3;break;
    case 3:  tvs2 = 4;break;
    case 4:  tvs2 = 3;break;
        }
    //    printf("in tvs2() ref==%d tvs2=%d\n", refBase, tvs2);
    return(tvs2);
  }
};

class JointGenoLk
{
public:
  double g11;
  double g12;
  double g22;
  double post11;
  double post12;
  double post22;
public: 
  JointGenoLk(){g11=g12=g22=post11=post12=post22=0.0;}
  void multiplyParentLikelihood(double lk);
  void CalcPost();
};

class JointGenoLk_denovo
{
public:
  double geno[10];
  double post[10];
public: 
  JointGenoLk_denovo(){
  for(int i=0; i<10; i++)
   {
    geno[i] = 1.0;
    post[i] = 0.0;
   }
  }
  void ResetGenoLk(double);
  void ResetPostLk(double);
  void multiplyParentLikelihood(double lk);
  void CalcPost();
};


//Def of pedigreeGLF class
//This class has the pedigree information and also the GLF handlers for each person
class PedigreeGLF
{
private:
  int nFam;
  int nPerson; 
  String glfFileKey;
  glfHandler * nonNULLglf;
  int nonNullIndex_i, nonNullIndex_j;  
  String nonNullPID;

 public:
  unsigned char refBase;
  int currentPos;
  int nFounders;
  Pedigree * ped;
  StringMap * glfMap;
  glfHandler **glf;
  std::vector<std::vector<int> > sexes;  
  int maleFounders;
  int femaleFounders;
 public:
  PedigreeGLF();
  PedigreeGLF(Pedigree *);
  ~PedigreeGLF();
  void InitializeGLFHandler(Pedigree *);
  glfHandler * GetNonNULLglf();
  void SetPedGLF(Pedigree * );
  void SetGLFMap(StringMap *);
  glfHandler ** GetPedGLFPointer();
  int GetFamCount();
  int GetFounderCout();
  int GetPersonCount();
  void GetSexes();
  bool Move2NextSection();
  bool Move2NextEntry();
  bool Move2NextBaseEntry();
  bool CheckSectionLabels();
  bool CheckSectionLabels(String&, int, int);
  //bool CheckConsistency();
  unsigned char GetRefBase();
};

#endif
