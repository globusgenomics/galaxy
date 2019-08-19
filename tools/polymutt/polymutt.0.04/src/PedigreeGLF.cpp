#include "PedigreeGLF.h"
#include <string.h>

#define MALE 1
#define FEMALE 2

void JointGenoLk::multiplyParentLikelihood(double lk)
{
  g11*=lk; g12*=lk; g22*=lk;
}

void JointGenoLk::CalcPost(){
  double sum = g11+g12+g22;
  if(sum==0.0) {
    post11=post12=post22=0.0; 
    return;
  }
  post11 = g11/sum;
  post12 = g12/sum;
  post22 = g22/sum;
}

void JointGenoLk_denovo::ResetGenoLk(double lk){
 for(int i=0; i<10; i++)
  geno[i] = lk;
}

void JointGenoLk_denovo::ResetPostLk(double lk){
 for(int i=0; i<10; i++)
  post[i] = lk;
}

void JointGenoLk_denovo::multiplyParentLikelihood(double lk)
{
  for(int i=0; i<10; i++)
   geno[i] *= lk;
}

void JointGenoLk_denovo::CalcPost(){
  double sum = 0.0;
  for(int i=0; i<10; i++)
   sum += geno[i];

  if(sum==0.0) {
    for(int i=0; i<10; i++)
     post[i] = 0.0; 
    return;
  }
  for(int i=0; i<10; i++)
    post[i] = geno[i]/sum;

}

PedigreeGLF::PedigreeGLF()
{
  nFam=0;
  nPerson = 0;
  nFounders = 0;
  currentPos = 0;
  glfMap = NULL;
  ped = NULL;
  glf = NULL;
  nonNULLglf = NULL;
  nonNullIndex_i = 0;
  nonNullIndex_j = 0;
  maleFounders=0;
  femaleFounders=0;
}

PedigreeGLF::PedigreeGLF(Pedigree * ped)
{
  SetPedGLF(ped);
}

PedigreeGLF::~PedigreeGLF()
{
  for(int i=0; i<nFam; i++)
    delete [] glf[i]; 
}

int PedigreeGLF::GetFamCount()
{
  return(nFam);
}
int PedigreeGLF::GetFounderCout()
{
  return(nFounders);
}
int PedigreeGLF::GetPersonCount()
{
 return(nPerson);
}

void PedigreeGLF::GetSexes()
{
 sexes.resize(nFam);
 for(int i=0; i<nFam; i++)
 {
  sexes[i].resize(ped->families[i]->count);
  for(int j=0; j<ped->families[i]->count; j++)
  {
   Person *p = ped->persons[ped->families[i]->path[j]];
   sexes[i][j] = p->sex;
   if(p->sex==MALE && p->isFounder()) maleFounders++;
   if(p->sex==FEMALE && p->isFounder()) femaleFounders++;
  }
 }
}

//prepare GLF files so that pointers of GLF are easily passed to functions for calculating likelihood
void PedigreeGLF::SetPedGLF(Pedigree * pedpt)
{
  int nValidGLF = 0;
  ped = pedpt;
  InitializeGLFHandler(ped);

  for(int i=0; i<ped->familyCount;i++) //iterate all families 
    {     
      nValidGLF = 0; 
      for(int j=0; j<ped->families[i]->count;j++) // iterate all members in a family
        {
          int idx = ped->families[i]->path[j];
          int glfIdx =  ped->persons[idx]->GetTraitID("GLF_Index");
          int glfFileIdx = (int)(ped->persons[idx]->traits[glfIdx]);
          glfFileKey.printf("%d", glfFileIdx);

	  if(glfFileIdx==0) { glf[i][j].handle=NULL;continue; }

          if(glfMap->Find(glfFileKey)<0)
          {
            warning("No entry found for the glf with the key [%s]\n", glfFileKey.c_str());
	  			glf[i][j].handle=NULL; continue;
	       }
          String glfFileName = * ((String*)((*glfMap).Object(glfFileKey)));

          bool flag = glf[i][j].Open(glfFileName.c_str()); 	 

	  if(nonNULLglf==NULL && flag) { 
		nonNULLglf = &glf[i][j];
		nonNullIndex_i = i;
		nonNullIndex_j = j;
		Person *p = ped->persons[ped->families[i]->path[j]];
		nonNullPID = p->pid;
	  }

          if(flag==false) 
	    error("GLF file %s can  not be opened!\n", glfFileName.c_str());
	  nValidGLF++;
        }
	if(nValidGLF==0)
		fprintf(stderr, "WARNING: No GLF files provided for family %s\n", ped->families[i]->famid.c_str());
    }
  maleFounders=0;
  femaleFounders=0;
  GetSexes();
}

void PedigreeGLF::SetGLFMap(StringMap * map)
{
  glfMap = map;
}

void PedigreeGLF::InitializeGLFHandler(Pedigree * ped)
{
  currentPos = 0;
  nFam = ped->familyCount;
  nPerson = ped->count;
  nonNULLglf = NULL;

  //Calculate the number of founders
  nFounders = 0;
  for(int i=0; i<ped->familyCount;i++) {
    int tmp = ped->families[i]->founders;
    nFounders += tmp;
  }
  
  if(glf!=NULL) delete [] glf;
  glf = new glfHandler * [nFam];
  for(int i=0; i<nFam; i++) {
    glf[i] = new glfHandler[ped->families[i]->count]; 
    if(glf[i]==NULL)
      error("glf mem alloc failed!\n"); 
  }
}

glfHandler * PedigreeGLF::GetNonNULLglf() { return(nonNULLglf); }

bool PedigreeGLF::Move2NextSection()
{
  bool flag = false; 
  for(int i=0; i<nFam; i++)
    {
      for(int j=0; j<ped->families[i]->count; j++)
        { 
	  if(glf[i][j].handle==NULL) continue;
	  if(glf[i][j].isStub) continue;
          flag = glf[i][j].NextSection(); 
	  if(glf[i][j].maxPosition!=nonNULLglf->maxPosition || glf[i][j].label !=nonNULLglf->label)
	   {
		Person *p = ped->persons[ped->families[i]->path[j]];
        	String pid = p->pid;
	   error("GLF files are not compatible:\n\tFile of person %s has section %s with %d entries ...\n\tFile of person %s has section %s with %d entries ...\n",
		nonNullPID.c_str(), glf[nonNullIndex_i][nonNullIndex_j].label.c_str(), 
		glf[nonNullIndex_i][nonNullIndex_j].maxPosition, pid.c_str(), glf[i][j].label.c_str(), glf[i][j].maxPosition);
	  }
          if(!flag) return(flag);
        }
    }
  currentPos = 0;
  return(flag);
}

bool PedigreeGLF::CheckSectionLabels()
{
  for(int i=0; i<nFam; i++)
    {
      for(int j=0; j<ped->families[i]->count; j++)
	if(glf[i][j].label.IsEmpty() || strchr(WHITESPACE,glf[i][j].label[0])!=NULL ) return(false);
    }
   return(true);
}

bool PedigreeGLF::CheckSectionLabels(String& pid, int f_idx, int p_idx)
{
  for(int i=0; i<nFam; i++)
    {
      for(int j=0; j<ped->families[i]->count; j++)
	if(glf[i][j].label.IsEmpty() || strchr(WHITESPACE,glf[i][j].label[0])!=NULL ){
	Person *p = ped->persons[ped->families[i]->path[j]];
	pid = p->pid;
	f_idx = i;
	p_idx = j;
	return(false);
	}
    }
   return(true);
}


bool PedigreeGLF::Move2NextEntry()
{
  bool flag = true;
  for(int i=0; i<nFam; i++)
    {
      for(int j=0; j<ped->families[i]->count; j++)
        {
	  if(glf[i][j].handle==NULL) continue;
          if(glf[i][j].position==currentPos){
            glf[i][j].NextEntry(); 
          }
        }
    }

  currentPos = nonNULLglf->position;
  refBase = nonNULLglf->data.refBase;
  for(int i=0; i<nFam; i++)
    {
      for(int j=0; j<ped->families[i]->count; j++)
        {
  	  if(glf[i][j].handle==NULL) continue;
          if(glf[i][j].position<currentPos){
            currentPos = glf[i][j].position;
            refBase = glf[i][j].data.refBase;
          }
        }
    }
  if(currentPos>nonNULLglf->maxPosition)
    flag = false;
  return(flag);
}


bool PedigreeGLF::Move2NextBaseEntry()
{
         if (currentPos > 0)
            {
            // Check whether we have reached the end of the current chromosome
            bool done = false;
  for(int i=0; i<nFam; i++)
      for(int j=0; j<ped->families[i]->count; j++)
	if (glf[i][j].handle!=NULL && glf[i][j].data.recordType == 0)
            return(false);
    }

  bool flag = true;
  for(int i=0; i<nFam; i++)
    {
      for(int j=0; j<ped->families[i]->count; j++)
        {
	  if(glf[i][j].handle==NULL) continue;
          if(glf[i][j].position==currentPos){
            glf[i][j].NextBaseEntry(); 
          }
        }
    }

  currentPos = nonNULLglf->position;
  refBase = nonNULLglf->data.refBase;

  for(int i=0; i<nFam; i++)
    {
      for(int j=0; j<ped->families[i]->count; j++)
        {
	  if(glf[i][j].handle==NULL) continue;
          if(glf[i][j].position<currentPos){
            currentPos = glf[i][j].position;
            refBase = glf[i][j].data.refBase;
          }
        }
    }
  if(currentPos>nonNULLglf->maxPosition)
    flag = false;

  return(flag);
}


unsigned char PedigreeGLF::GetRefBase()
{
  return(refBase);
}

glfHandler ** PedigreeGLF::GetPedGLFPointer()
{
  return(glf);
}
