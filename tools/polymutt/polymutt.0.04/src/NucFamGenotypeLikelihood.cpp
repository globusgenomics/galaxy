#include "NucFamGenotypeLikelihood.h"
#include <vector>

void NucFamGenotypeLikelihood::InitializeValues()
{
  vcf = false;
  nFam = 0;
  nPerson = 0;
  nFounders = 0;
  prior = 0.0;
  priorFreq = 0.0;
  theta = 0.001;
  totalDepth = 0;
  avgMapQual = 0.0;
  avgDepth = 0.0;
  numSampWithData=0;
  percSampWithData = 0.0;
  currentFam = 0;
  usePriorFreq = true;
  useMleFreq = false;
  pedGLF = NULL;
  kidsCondLikelihood = NULL;
  transMissionProb = NULL;
  parentMarginal = NULL;
  parentConditional = NULL;
  parentGLF = NULL;
  postProb = NULL ;
  dosage = NULL;
  bestGenoIdx = NULL;
  bestGenoLabel = NULL;  
  varllk.resize(8);
  varllk_noprior.resize(8);
  varfreq.resize(8);
  denovo_mono = false;
  denovoLR = -1;
}

NucFamGenotypeLikelihood::NucFamGenotypeLikelihood()
{
  InitializeValues();

}

NucFamGenotypeLikelihood::NucFamGenotypeLikelihood(PedigreeGLF * pedglf)
{
  InitializeValues();
  SetGLF(pedglf);
}

NucFamGenotypeLikelihood::~NucFamGenotypeLikelihood()
{
  DeleteTransmissionProb();
  DeletePostProb();
  DeleteParentMarginal();  
}

void NucFamGenotypeLikelihood::SetGLF(PedigreeGLF *pedglf)
{
  pedGLF = pedglf;
  ped = pedglf->ped;
  nFam = pedglf->GetFamCount();
  nPerson = pedGLF->ped->count;
  nFounders = pedGLF->nFounders;
  InitializeTransmissionProb();
  InitializeParentMarginal();
  InitializePostProb();
  SetPolyPrior();
}

void NucFamGenotypeLikelihood::SetPolyQual(double qual)
{ polyQual = qual;}

void NucFamGenotypeLikelihood::SetMleFreq(double freq)
{ min = freq; }

void  NucFamGenotypeLikelihood::SetAlleles(int a1, int a2) 
{ 
  allele1 = a1; 
  allele2 = a2;
  geno11 = glfHandler::GenotypeIndex(a1, a1);
  geno12 = glfHandler::GenotypeIndex(a1, a2);
  geno22 = glfHandler::GenotypeIndex(a2, a2);

}

void NucFamGenotypeLikelihood::SetDenovoMutationModel()
{
  aM.SetMutationRate(par->denovo_mut_rate);
  aM.SetTsTvRatio(par->denovo_tstv_ratio);
  aM.SetAlleleMutMatrix();
  gM.SetGenoMutMatrix(aM);
}

void NucFamGenotypeLikelihood::InitializeTransmissionProb()
{
  transMissionProb = new double *[9];
  for(int i=0; i<9; i++) { 
    transMissionProb[i] = new double[3]; 
    for(int j=0; j<3; j++)
      transMissionProb[i][j] = 0.0;
  }
  transMissionProb[0][0] = 1.0;
  transMissionProb[1][0] = 0.5;  transMissionProb[1][1] = 0.5;
  transMissionProb[2][1] = 1.0;
  transMissionProb[3][0] = 0.5;  transMissionProb[3][1] = 0.5;
  transMissionProb[4][0] = 0.25; transMissionProb[4][1] = 0.5; transMissionProb[4][2] = 0.25;
  transMissionProb[5][1] = 0.5;  transMissionProb[5][2] = 0.5;
  transMissionProb[6][1] = 1.0;
  transMissionProb[7][1] = 0.5;  transMissionProb[7][2] = 0.5;
  transMissionProb[8][2] = 1.0;
}

void NucFamGenotypeLikelihood::DeleteTransmissionProb()
{
  for(int i=0; i<9; i++) { 
    if(transMissionProb[i]!=NULL)
      delete [] transMissionProb[i];
  }
}

void NucFamGenotypeLikelihood::InitializeKidsCondLikelihood()
{
  kidsCondLikelihood = new double **[nFam];
  for(int i=0 ; i<nFam; i++){
    int famSize = pedGLF->ped->families[i]->count;
    kidsCondLikelihood[i] = new double *[famSize];
    for(int j=0; j<pedGLF->ped->families[i]->count; j++){
      kidsCondLikelihood[i][j] = new double [9];
      for(int k=0; k<9; k++)
	kidsCondLikelihood[i][j][k] = 0.0;
    }
  }
}


void NucFamGenotypeLikelihood::DeleteKidsCondLikelihood()
{
  for(int i=0 ; i<nFam; i++)
    {   
      for(int j=0; j<pedGLF->ped->families[i]->count; j++)
	delete [] kidsCondLikelihood[i][j];
      delete kidsCondLikelihood[i];
    }
}

void NucFamGenotypeLikelihood::InitializeParentMarginal()
{
  parentMarginal = new double *[nFam];
  parentConditional = new double *[nFam];
  parentGLF      = new double *[nFam];
  for(int i=0; i<nFam; i++)
    {
      parentMarginal[i] = new double[9];
      parentConditional[i] = new double[9];
      parentGLF[i] = new double[9];
      for(int j=0; j<9; j++)
	{
	  parentMarginal[i][j] = 0.0;
          parentConditional[i][j] = 0.0;
	  parentGLF[i][j] = 0.0;
	}
    }
}

void NucFamGenotypeLikelihood::DeleteParentMarginal()
{
  for(int i=0; i<nFam; i++) {
    delete [] parentMarginal[i];
    delete [] parentConditional[i];
    delete [] parentGLF[i];
  }
}

void NucFamGenotypeLikelihood::InitializePostProb()
{
  postProb = new double **[nFam];
  bestGenoIdx = new int *[nFam];
  bestGenoLabel = new String *[nFam];
  dosage = new double *[nFam];
  for(int i=0; i<nFam; i++)
    {
      postProb[i] = new double *[pedGLF->ped->families[i]->count];
      bestGenoIdx[i] = new int[pedGLF->ped->families[i]->count];
      bestGenoLabel[i] = new String[pedGLF->ped->families[i]->count];
      dosage[i] = new double[pedGLF->ped->families[i]->count];

      for(int j=0; j<pedGLF->ped->families[i]->count; j++)
	{
	  postProb[i][j] = new double[10];
	  bestGenoIdx[i][j] = 0;
	  bestGenoLabel[i][j] = "";
	  for(int k=0; k<10; k++)
	    postProb[i][j][k] = 0.0;
	}
    }
}

void NucFamGenotypeLikelihood::DeletePostProb()
{
  for(int i=0; i<nFam; i++){
    for(int j=0; j<pedGLF->ped->families[i]->count; j++){
      delete [] postProb[i][j];
    }
    delete [] bestGenoIdx[i];
    delete [] bestGenoLabel[i];
    delete [] postProb[i];
  }
}
void NucFamGenotypeLikelihood::SetTheta(double thetaValue)
{theta = thetaValue;}

void NucFamGenotypeLikelihood::SetPolyPrior()
{
  prior = 0;
  if(nFounders==0)
    error("Family size is zero\n");
  
  for(int i=1; i<=2*nFounders; i++) //number of chromosomes
    prior += 1.0 / i;
  priorFreq = 1-1./prior;
  prior *= theta;  
  
}

double NucFamGenotypeLikelihood::GetPolyPrior() { return(prior); }
double NucFamGenotypeLikelihood::GetPriorFreq() { return priorFreq; }
void NucFamGenotypeLikelihood::SetMleFreqFlag(bool b) { useMleFreq = b; }
void NucFamGenotypeLikelihood::SetPriorFreqFlag(bool b) { usePriorFreq = b; }

void NucFamGenotypeLikelihood::SetParentPrior(double freq)
{
  if(nFam>1) {
    parentPrior[0] = pow(freq,4);
    parentPrior[1] = freq*freq * freq*(1-freq)*2;
    parentPrior[2] = freq*freq * (1-freq)*(1-freq);
    parentPrior[3] = freq*(1-freq)*2 * freq*freq;
    parentPrior[4] = freq*(1-freq)*2 * freq*(1-freq)*2;
    parentPrior[5] = freq*(1-freq)*2 * (1-freq)*(1-freq);
    parentPrior[6] = (1-freq)*(1-freq) * freq*freq;
    parentPrior[7] = (1-freq)*(1-freq) * freq*(1-freq)*2;
    parentPrior[8] = (1-freq)*(1-freq) * (1-freq)*(1-freq);
  }
  else SetParentPriorSingleTrio();
}

void NucFamGenotypeLikelihood::SetParentPrior_denovo(double freq)
{
    parentPrior[0] = pow(freq,4);
    parentPrior[1] = freq*freq * freq*(1-freq)*2;
    parentPrior[2] = freq*freq * (1-freq)*(1-freq);
    parentPrior[3] = freq*(1-freq)*2 * freq*freq;
    parentPrior[4] = freq*(1-freq)*2 * freq*(1-freq)*2;
    parentPrior[5] = freq*(1-freq)*2 * (1-freq)*(1-freq);
    parentPrior[6] = (1-freq)*(1-freq) * freq*freq;
    parentPrior[7] = (1-freq)*(1-freq) * freq*(1-freq)*2;
    parentPrior[8] = (1-freq)*(1-freq) * (1-freq)*(1-freq);
}

void NucFamGenotypeLikelihood::SetParentPriorSingleTrio()
{
  parentPrior[0] = 0.0;
  parentPrior[1] = 0.24;
  parentPrior[2] = 0.04;
  parentPrior[3] = 0.24;
  parentPrior[4] = 0.16;
  parentPrior[5] = 0.08;
  parentPrior[6] = 0.04;
  parentPrior[7] = 0.08;
  parentPrior[8] = 0.12;
}

void NucFamGenotypeLikelihood::SetParentPriorSingleTrio_denovo(double freq)
{
  if(freq!=1.0) {
  parentPrior[0] = 0.0;
  parentPrior[1] = 0.24;
  parentPrior[2] = 0.04;
  parentPrior[3] = 0.24;
  parentPrior[4] = 0.16;
  parentPrior[5] = 0.08;
  parentPrior[6] = 0.04;
  parentPrior[7] = 0.08;
  parentPrior[8] = 0.12;
  }
  else {
    parentPrior[0] = pow(freq,4);
    parentPrior[1] = freq*freq * freq*(1-freq)*2;
    parentPrior[2] = freq*freq * (1-freq)*(1-freq);
    parentPrior[3] = freq*(1-freq)*2 * freq*freq;
    parentPrior[4] = freq*(1-freq)*2 * freq*(1-freq)*2;
    parentPrior[5] = freq*(1-freq)*2 * (1-freq)*(1-freq);
    parentPrior[6] = (1-freq)*(1-freq) * freq*freq;
    parentPrior[7] = (1-freq)*(1-freq) * freq*(1-freq)*2;
    parentPrior[8] = (1-freq)*(1-freq) * (1-freq)*(1-freq);
  }
}

int NucFamGenotypeLikelihood::GetFamCount()
{ return(nFam); }

double NucFamGenotypeLikelihood::f(double freq)
{
  if(par->denovo)
    return(-(allFamLogLikelihood_denovo(freq)));
  return(-(allFamLogLikelihood(freq)));
}

double NucFamGenotypeLikelihood::OptimizeFrequency()
{
  a = 0.0001; 
  fa = f(a);
  b = 0.9999; 
  fb = f(b);
  c = 0.5; fc = f(c);

  Brent(par->precision);

  return min;
}

double NucFamGenotypeLikelihood::GetMaxLogLikelihood(){return(-fmin);}
double NucFamGenotypeLikelihood::GetMaxLikelihood(){return(pow10(-fmin));}
double NucFamGenotypeLikelihood::GetMinimizer(){return(min);}

double NucFamGenotypeLikelihood::PolymorphismLikelihood(int  a1, int a2)
{
  SetAlleles(a1, a2);
  
  if(nFam>1)
    {
      OptimizeFrequency();
      return(GetMaxLikelihood());
    }
  else {
    double dummy = 1.0;
    return(allFamLikelihood(dummy));
  }
}


double NucFamGenotypeLikelihood::MonomorphismLikelihood(int refBase)
{
  double lRef = 1.0;
  const double *ptrlk;
  
  int homoRefIdx = glfHandler::GenotypeIndex(pedGLF->refBase, pedGLF->refBase);
  for(int i=0; i<nFam; i++)
    for(int j=0; j<pedGLF->ped->families[i]->count; j++)
      {
	if(pedGLF->glf[i][j].handle==NULL)continue;
	ptrlk = pedGLF->glf[i][j].GetLikelihoods(pedGLF->currentPos);
	lRef *= ptrlk[homoRefIdx]; 
      }
  return(lRef);
}

double NucFamGenotypeLikelihood::PolymorphismLogLikelihood(int  a1, int a2)
{
  SetAlleles(a1, a2);
  if(nFam>1)
    {
      OptimizeFrequency();
      return(GetMaxLogLikelihood());
    }
  else {
    double dummy = 1.0;
    return(allFamLogLikelihood(dummy));
  }
}

double NucFamGenotypeLikelihood::MonomorphismLogLikelihood_denovo(int  refBase, int a2)
{
  SetAlleles(refBase, a2);
  return(allFamLogLikelihood(1.0));
}

double NucFamGenotypeLikelihood::MonomorphismLogLikelihood(int refBase)
{
  double lRef = 0.0;
  const unsigned char * ptrlk;
  int i, j;
  int homoRefIdx = glfHandler::GenotypeIndex(pedGLF->refBase, pedGLF->refBase);

  for(i=0; i<nFam; i++)
    for(j=0; j<pedGLF->ped->families[i]->count; j++)
      {
	if(pedGLF->glf[i][j].handle==NULL) continue;
	ptrlk = pedGLF->glf[i][j].GetLogLikelihoods(pedGLF->currentPos);
	lRef += -double(ptrlk[homoRefIdx])/10;
      }
  return(lRef);
}


void NucFamGenotypeLikelihood::CalcReadStats()
{
  totalDepth = 0;
  numSampWithData = 0;
  avgMapQual = 0.0;
  avgDepth = 0.0;
  
  for(int i=0; i<nFam; i++)
    for(int j=0; j<pedGLF->ped->families[i]->count; j++)
      {
	if(pedGLF->glf[i][j].handle==NULL) continue;
	totalDepth +=  pedGLF->glf[i][j].GetDepth(pedGLF->currentPos); 
	avgMapQual += pedGLF->glf[i][j].GetMapQuality(pedGLF->currentPos);
	if(pedGLF->glf[i][j].GetDepth(pedGLF->currentPos)>0) numSampWithData++;
      }
  
  if(numSampWithData==0) {
    avgDepth = 0.;
    avgMapQual = 0.;
    percSampWithData = 0.;
  }
  else {
    avgDepth = double(totalDepth)/double(numSampWithData);
    avgMapQual /= double(numSampWithData);
    percSampWithData = (double(numSampWithData)/double(nPerson));
  }
}

void NucFamGenotypeLikelihood::CalcDosage(int a1, int a2, double freq)
{
  CalcPostProb(freq);
}

void NucFamGenotypeLikelihood::GetBestGenotype(int a1, int a2, double freq)
{}

void NucFamGenotypeLikelihood::CalcPostProb(double freq)
{
  for(int i=0; i<nFam; i++)
   if(par->denovo==false)
    CalcPostProb_SingleNucFam(i, freq);
   else
    CalcPostProb_SingleNucFam_denovo(i, freq);
}

void NucFamGenotypeLikelihood::CalcPostProb_SingleNucFam(int i, double freq)
{
  double p11, p12, p22;
  double sum;
  int best;

  if(pedGLF->ped->families[i]->count<=pedGLF->ped->families[i]->founders)
  {
   for(int j=0; j<pedGLF->ped->families[i]->founders; j++)
    CalcPostProb_SinglePerson(i, j, freq);
    return;
  }

      CalcParentMarginal(i, freq);
      
      for(int j=0; j<pedGLF->ped->families[i]->count; j++)
	   {
	    if(j==0){
	    p11 = parentMarginal[i][0] + parentMarginal[i][1] + parentMarginal[i][2];
	    p12 = parentMarginal[i][3] + parentMarginal[i][4] + parentMarginal[i][5];
	    p22 = parentMarginal[i][6] + parentMarginal[i][7] + parentMarginal[i][8];
	    sum = p11 + p12 + p22;

	    if(sum==0) 
		   postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 1/3;
	    else {
	      postProb[i][j][0] = p11/sum;
	      postProb[i][j][1] = p12/sum;
	      postProb[i][j][2] = p22/sum;
	    }
	    best = GetBestGenoIdx(p11, p12, p22);
	    bestGenoIdx[i][j] = best; 
	    bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best); 
	    //bestGenoLabel[i][j] = GetBestGenoLabel(best); 
	    dosage[i][j] = CalcDosage(i,j); 
	  } //first parent
	  else if(j==1){
	    p11 = parentMarginal[i][0] + parentMarginal[i][3] + parentMarginal[i][6];
	    p12 = parentMarginal[i][1] + parentMarginal[i][4] + parentMarginal[i][7];
	    p22 = parentMarginal[i][2] + parentMarginal[i][5] + parentMarginal[i][8];
	    sum = p11 + p12 + p22;
            if(sum==0) 
                postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 1/3;
	   else {
	      postProb[i][j][0] = p11/sum;
	      postProb[i][j][1] = p12/sum;
	      postProb[i][j][2] = p22/sum;
	    }
	    best = GetBestGenoIdx(p11, p12, p22);
	    bestGenoIdx[i][j] = best; 
	    bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best);
	    //bestGenoLabel[i][j] = GetBestGenoLabel(best);
	    dosage[i][j] = CalcDosage(i,j);
	  } //second parent
	  else {

	    JointGenoLk JGLK = KidJointGenoLikelihood(i,j); //Joint likelihood of the kid j marginalizing other kids
	    
	    JGLK.CalcPost();
	    
	    postProb[i][j][0] = JGLK.post11;
	    postProb[i][j][1] = JGLK.post12;
	    postProb[i][j][2] = JGLK.post22;
	    
	    best = GetBestGenoIdx(JGLK.post11, JGLK.post12, JGLK.post22);

	    bestGenoIdx[i][j] = best;
	    bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best);
	    //bestGenoLabel[i][j] = GetBestGenoLabel(best);
	    dosage[i][j] = CalcDosage(i,j);

	  }  //kids
      } //END of for loop
}

void NucFamGenotypeLikelihood::CalcPostProb_SingleNucFam_denovo(int i, double freq)
{
  double p11, p12, p22;
  double sum;
  int best;

  if(pedGLF->ped->families[i]->count<=pedGLF->ped->families[i]->founders)
    {
      for(int j=0; j<pedGLF->ped->families[i]->founders; j++)
	CalcPostProb_SinglePerson(i, j, freq);
      return;
    }
  
  CalcParentMarginal_denovo(i, freq);
  
  for(int j=0; j<pedGLF->ped->families[i]->count; j++)
    {
      if(j==0){
	p11 = parentMarginal[i][0] + parentMarginal[i][1] + parentMarginal[i][2];
	p12 = parentMarginal[i][3] + parentMarginal[i][4] + parentMarginal[i][5];
	p22 = parentMarginal[i][6] + parentMarginal[i][7] + parentMarginal[i][8];
	sum = p11 + p12 + p22;
	
	if(sum==0) 
	  postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 1/3;
	else {
	  postProb[i][j][0] = p11/sum;
	  postProb[i][j][1] = p12/sum;
	  postProb[i][j][2] = p22/sum;
	}
	best = GetBestGenoIdx(p11, p12, p22);
	bestGenoIdx[i][j] = best; 
	//bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best); 
	bestGenoLabel[i][j] = GetBestGenoLabel(best); 
	dosage[i][j] = CalcDosage(i,j); 
      } //first parent
      else if(j==1){
	p11 = parentMarginal[i][0] + parentMarginal[i][3] + parentMarginal[i][6];
	p12 = parentMarginal[i][1] + parentMarginal[i][4] + parentMarginal[i][7];
	p22 = parentMarginal[i][2] + parentMarginal[i][5] + parentMarginal[i][8];
	sum = p11 + p12 + p22;
	if(sum==0) 
	  postProb[i][j][0] = postProb[i][j][1] = postProb[i][j][2] = 1/3;
	else {
	  postProb[i][j][0] = p11/sum;
	  postProb[i][j][1] = p12/sum;
	  postProb[i][j][2] = p22/sum;
	}
	best = GetBestGenoIdx(p11, p12, p22);
	bestGenoIdx[i][j] = best; 
	//bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best);
	bestGenoLabel[i][j] = GetBestGenoLabel(best);
	dosage[i][j] = CalcDosage(i,j);
      } //second parent
      else {
	
	JointGenoLk_denovo JGLK = KidJointGenoLikelihood_denovo(i,j); //Joint genotype likelihood of the kid j marginalizing other kids
	
	JGLK.CalcPost();
	
	for(int k=0; k<10; k++)
	  postProb[i][j][k] = JGLK.post[k];
	
	double maxPost = 0.0;
	best = 0;
	for(int k=0; k<10; k++)
	  {
	    if(maxPost<JGLK.post[k])
	      {
		maxPost = JGLK.post[k];
		best = k;
	      }
	  }
	
	bestGenoIdx[i][j] = best;
	//bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best);
	bestGenoLabel[i][j] = GetBestGenoLabel_denovo(best);
	//dosage[i][j] = CalcDosage(i,j);
	dosage[i][j] = 0.0;
      }  //kids
    } //END of for loop
}

void NucFamGenotypeLikelihood::CalcPostProb_SinglePerson(int i, int j, double freq)
{
  double lk11, lk12, lk22;
  getGenoLikelihood(&(pedGLF->glf[i][j]), allele1, allele2, &lk11, &lk12, &lk22);
  double mlk11 = lk11*freq*freq;
  double mlk12 = lk12*freq*(1-freq)*2;
  double mlk22 = lk22*(1-freq)*(1-freq);
  double sum = mlk11+mlk12+mlk22;
  if(sum==0) postProb[i][j][0]=postProb[i][j][1]=postProb[i][j][2]=1/3;
  else {
    postProb[i][j][0] =  mlk11/sum;
    postProb[i][j][1] =  mlk12/sum;
    postProb[i][j][2] =  mlk22/sum;
  }
  
  int best = GetBestGenoIdx(mlk11, mlk12, mlk22);
  bestGenoIdx[i][j] =  best;
  bestGenoLabel[i][j] = GetBestGenoLabel_vcfv4(best);
  //bestGenoLabel[i][j] = GetBestGenoLabel(best);
  dosage[i][j] = CalcDosage(i, j);
}

// Likelihood of genotypes of a kid marginalizing parents and (s)his siblings
JointGenoLk NucFamGenotypeLikelihood::KidJointGenoLikelihood(int famIdx, int kidIdx)
{
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22;

  getGenoLikelihood(pedGLF->glf[famIdx], allele1, allele2, &lkF11, &lkF12, &lkF22);
  getGenoLikelihood(pedGLF->glf[famIdx]+1, allele1, allele2, &lkM11, &lkM12, &lkM22);
  
  JointGenoLk JGLK1111 = likelihoodKidGenotype(allele1,allele1,allele1,allele1, famIdx, kidIdx);
  JointGenoLk JGLK1112 = likelihoodKidGenotype(allele1,allele1,allele1,allele2, famIdx, kidIdx);
  JointGenoLk JGLK1122 = likelihoodKidGenotype(allele1,allele1,allele2,allele2, famIdx, kidIdx);
  JointGenoLk JGLK1211 = likelihoodKidGenotype(allele1,allele2,allele1,allele1, famIdx, kidIdx);
  JointGenoLk JGLK1212 = likelihoodKidGenotype(allele1,allele2,allele1,allele2, famIdx, kidIdx);
  JointGenoLk JGLK1222 = likelihoodKidGenotype(allele1,allele2,allele2,allele2, famIdx, kidIdx);
  JointGenoLk JGLK2211 = likelihoodKidGenotype(allele2,allele2,allele1,allele1, famIdx, kidIdx);
  JointGenoLk JGLK2212 = likelihoodKidGenotype(allele2,allele2,allele1,allele2, famIdx, kidIdx);
  JointGenoLk JGLK2222 = likelihoodKidGenotype(allele2,allele2,allele2,allele2, famIdx, kidIdx);
  
  JGLK1111.multiplyParentLikelihood(parentGLF[famIdx][0]*parentPrior[0]);
  JGLK1112.multiplyParentLikelihood(parentGLF[famIdx][1]*parentPrior[1]);
  JGLK1122.multiplyParentLikelihood(parentGLF[famIdx][2]*parentPrior[2]);
  JGLK1211.multiplyParentLikelihood(parentGLF[famIdx][3]*parentPrior[3]);
  JGLK1212.multiplyParentLikelihood(parentGLF[famIdx][4]*parentPrior[4]);
  JGLK1222.multiplyParentLikelihood(parentGLF[famIdx][5]*parentPrior[5]);
  JGLK2211.multiplyParentLikelihood(parentGLF[famIdx][6]*parentPrior[6]);
  JGLK2212.multiplyParentLikelihood(parentGLF[famIdx][7]*parentPrior[7]);
  JGLK2222.multiplyParentLikelihood(parentGLF[famIdx][8]*parentPrior[8]);
  
  double JGLK11 = JGLK1111.g11 + JGLK1112.g11 + JGLK1122.g11 + JGLK1211.g11 + JGLK1212.g11 + JGLK1222.g11 + JGLK2211.g11 + JGLK2212.g11 + JGLK2222.g11;
  double JGLK12 = JGLK1111.g12 + JGLK1112.g12 + JGLK1122.g12 + JGLK1211.g12 + JGLK1212.g12 + JGLK1222.g12 + JGLK2211.g12 + JGLK2212.g12 + JGLK2222.g12;
  double JGLK22 = JGLK1111.g22 + JGLK1112.g22 + JGLK1122.g22 + JGLK1211.g22 + JGLK1212.g22 + JGLK1222.g22 + JGLK2211.g22 + JGLK2212.g22 + JGLK2222.g22;
  
  JointGenoLk jglk;
  jglk.g11 = JGLK11;
  jglk.g12 = JGLK12;
  jglk.g22 = JGLK22;

  return(jglk);
}

// Likelihood of genotypes of a kid marginalizing parents and (s)his siblings allowing for de novo mutations
JointGenoLk_denovo NucFamGenotypeLikelihood::KidJointGenoLikelihood_denovo(int famIdx, int kidIdx)
{
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22;
  
  getGenoLikelihood(pedGLF->glf[famIdx], allele1, allele2, &lkF11, &lkF12, &lkF22);
  getGenoLikelihood(pedGLF->glf[famIdx]+1, allele1, allele2, &lkM11, &lkM12, &lkM22);
  
  JointGenoLk_denovo JGLK[9];
  JGLK[0] = likelihoodKidGenotype_denovo(allele1,allele1,allele1,allele1, famIdx, kidIdx);
  JGLK[1] = likelihoodKidGenotype_denovo(allele1,allele1,allele1,allele2, famIdx, kidIdx);
  JGLK[2] = likelihoodKidGenotype_denovo(allele1,allele1,allele2,allele2, famIdx, kidIdx);
  JGLK[3] = likelihoodKidGenotype_denovo(allele1,allele2,allele1,allele1, famIdx, kidIdx);
  JGLK[4] = likelihoodKidGenotype_denovo(allele1,allele2,allele1,allele2, famIdx, kidIdx);
  JGLK[5] = likelihoodKidGenotype_denovo(allele1,allele2,allele2,allele2, famIdx, kidIdx);
  JGLK[6] = likelihoodKidGenotype_denovo(allele2,allele2,allele1,allele1, famIdx, kidIdx);
  JGLK[7] = likelihoodKidGenotype_denovo(allele2,allele2,allele1,allele2, famIdx, kidIdx);
  JGLK[8] = likelihoodKidGenotype_denovo(allele2,allele2,allele2,allele2, famIdx, kidIdx);
  
  for(int i=0; i<9; i++)
    JGLK[i].multiplyParentLikelihood(parentGLF[famIdx][i]*parentPrior[i]);
  
  JointGenoLk_denovo jglk;  
  jglk.ResetGenoLk(0.0);
  
  for(int i=0; i<10; i++)
   for(int j=0; j<9; j++)
   {
    jglk.geno[i] += JGLK[j].geno[i];
   }
  return(jglk);
}

void NucFamGenotypeLikelihood::toOriginalLikelihood(int loglk11, int loglk12, int loglk22, double *lk11, double *lk12, double *lk22)
{
  *lk11 = pow10(-double(loglk11)/10);
  *lk12 = pow10(-double(loglk12)/10);
  *lk22 = pow10(-double(loglk22)/10);
}

double NucFamGenotypeLikelihood::toOriginalLikelihood(int loglk)
{
  double lk = pow10(-double(loglk)/10);
  return(lk);
}

double NucFamGenotypeLikelihood::allFamLikelihood(double freq)
{
  double loglk = 0.0;
  
# ifdef _OPENMP
# pragma omp parallel for reduction(+:loglk)
# endif
  for(int i=0; i<nFam; i++)
    {
      loglk += log10(lkSingleFam(i, freq));
    }
  return(pow10(loglk));
}

double NucFamGenotypeLikelihood::allFamLikelihood_denovo(double freq)
{
  double loglk = 0.0;
  
# ifdef _OPENMP
# pragma omp parallel for reduction(+:loglk)
# endif
  for(int i=0; i<nFam; i++)
    {
      loglk += log10(lkSingleFam_denovo(i, freq));
    }
  return(pow10(loglk));
}

double NucFamGenotypeLikelihood::allFamLogLikelihood(double freq)
{
  double loglk = 0.0;
  
# ifdef _OPENMP
# pragma omp parallel for reduction(+:loglk)
# endif
  for(int i=0; i<nFam; i++)
    {
      loglk += log10(lkSingleFam(i, freq));
    }
  
  return(loglk);
}

double NucFamGenotypeLikelihood::allFamLogLikelihood_denovo(double freq)
{
  double loglk = 0.0;
  
# ifdef _OPENMP
# pragma omp parallel for reduction(+:loglk)
# endif
  for(int i=0; i<nFam; i++)
    {
      loglk += log10(lkSingleFam_denovo(i, freq));
    }
  
  return(loglk);
}

double NucFamGenotypeLikelihood::lkSingleFam(int i, double freq)
{
  if(pedGLF->ped->families[i]->count==pedGLF->ped->families[i]->founders)
    {
   double lk = 1.0;
  	for(int j=0; j<pedGLF->ped->families[i]->founders; j++)
     lk *= lkSinglePerson(i, j, freq);
   return lk;
  }
  double sum = 0.0;
  
  CalcParentMarginal(i, freq);
  for(int idx=0; idx<9; idx++)
    sum += parentMarginal[i][idx]; //Note that Marginals here are joint likelihood of parents marginalizing offspring
  
  return(sum);
}

double NucFamGenotypeLikelihood::lkSingleFam_denovo(int i, double freq)
{
  if(pedGLF->ped->families[i]->count==pedGLF->ped->families[i]->founders)
    {
   double lk = 1.0;
  	for(int j=0; j<pedGLF->ped->families[i]->founders; j++)
     lk *= lkSinglePerson(i, j, freq);
   return lk;
  }
  double sum = 0.0;
  
  CalcParentMarginal_denovo(i, freq);
  for(int idx=0; idx<9; idx++)
    sum += parentMarginal[i][idx]; //Note that Marginals here are joint likelihood of parents marginalizing offspring
  
  return(sum);
}

double NucFamGenotypeLikelihood::logLkSingleFam(int i, double freq)
{
  return(log10(lkSingleFam(i, freq)));
}

double NucFamGenotypeLikelihood::logLkSingleFam_denovo(int i, double freq)
{
  return(log10(lkSingleFam_denovo(i, freq)));
}

double NucFamGenotypeLikelihood::lkSinglePerson(int i, int j, double freq)
{
  double sum=0.0;
  double lk11, lk12, lk22;
  getGenoLikelihood(&(pedGLF->glf[i][j]), allele1, allele2, &lk11, &lk12, &lk22);
  sum = sum + lk11*freq*freq + lk12*freq*(1-freq)*2 + lk22*(1-freq)*(1-freq);
  return(sum);
}

void NucFamGenotypeLikelihood::CalcParentMarginal(int i, double freq)
{
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22; // lkC11, lkC12, lkC22;
  lkF11=lkF12=lkF22=lkM11=lkM12=lkM22=1.0;

  getGenoLikelihood(pedGLF->glf[i], allele1, allele2, &lkF11, &lkF12, &lkF22);
  getGenoLikelihood(pedGLF->glf[i]+1, allele1, allele2, &lkM11, &lkM12, &lkM22);

  parentGLF[i][0] = lkF11*lkM11;
  parentGLF[i][1] = lkF11*lkM12;
  parentGLF[i][2] = lkF11*lkM22;
  parentGLF[i][3] = lkF12*lkM11;
  parentGLF[i][4] = lkF12*lkM12;
  parentGLF[i][5] = lkF12*lkM22;
  parentGLF[i][6] = lkF22*lkM11;
  parentGLF[i][7] = lkF22*lkM12;
  parentGLF[i][8] = lkF22*lkM22;


  if(nFam>1)
    SetParentPrior(freq);
  else 
    SetParentPriorSingleTrio();

  parentConditional[i][0] = likelihoodKids(allele1,allele1,allele1,allele1, i)*parentGLF[i][0];
  parentConditional[i][1] = likelihoodKids(allele1,allele1,allele1,allele2, i)*parentGLF[i][1];
  parentConditional[i][2] = likelihoodKids(allele1,allele1,allele2,allele2, i)*parentGLF[i][2];
  parentConditional[i][3] = likelihoodKids(allele1,allele2,allele1,allele1, i)*parentGLF[i][3];
  parentConditional[i][4] = likelihoodKids(allele1,allele2,allele1,allele2, i)*parentGLF[i][4];
  parentConditional[i][5] = likelihoodKids(allele1,allele2,allele2,allele2, i)*parentGLF[i][5];
  parentConditional[i][6] = likelihoodKids(allele2,allele2,allele1,allele1, i)*parentGLF[i][6];
  parentConditional[i][7] = likelihoodKids(allele2,allele2,allele1,allele2, i)*parentGLF[i][7];
  parentConditional[i][8] = likelihoodKids(allele2,allele2,allele2,allele2, i)*parentGLF[i][8];

  //Mariginal is joint likelihood of parents marginalizing all kids
  //This is to avoid redundant calculation
  parentMarginal[i][0] = parentConditional[i][0]*parentPrior[0];
  parentMarginal[i][1] = parentConditional[i][1]*parentPrior[1];
  parentMarginal[i][2] = parentConditional[i][2]*parentPrior[2];
  parentMarginal[i][3] = parentConditional[i][3]*parentPrior[3];
  parentMarginal[i][4] = parentConditional[i][4]*parentPrior[4];
  parentMarginal[i][5] = parentConditional[i][5]*parentPrior[5];
  parentMarginal[i][6] = parentConditional[i][6]*parentPrior[6];
  parentMarginal[i][7] = parentConditional[i][7]*parentPrior[7];
  parentMarginal[i][8] = parentConditional[i][8]*parentPrior[8];

}

void NucFamGenotypeLikelihood::CalcParentMarginal_denovo(int i, double freq)
{
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22; // lkC11, lkC12, lkC22;
  lkF11=lkF12=lkF22=lkM11=lkM12=lkM22=1.0;

  getGenoLikelihood(pedGLF->glf[i], allele1, allele2, &lkF11, &lkF12, &lkF22);
  getGenoLikelihood(pedGLF->glf[i]+1, allele1, allele2, &lkM11, &lkM12, &lkM22);

  parentGLF[i][0] = lkF11*lkM11;
  parentGLF[i][1] = lkF11*lkM12;
  parentGLF[i][2] = lkF11*lkM22;
  parentGLF[i][3] = lkF12*lkM11;
  parentGLF[i][4] = lkF12*lkM12;
  parentGLF[i][5] = lkF12*lkM22;
  parentGLF[i][6] = lkF22*lkM11;
  parentGLF[i][7] = lkF22*lkM12;
  parentGLF[i][8] = lkF22*lkM22;


  if(nFam>1)
    SetParentPrior_denovo(freq);
  else 
    SetParentPriorSingleTrio_denovo(freq);

  parentConditional[i][0] = likelihoodKids_denovo(allele1,allele1,allele1,allele1, i)*parentGLF[i][0];
  parentConditional[i][1] = likelihoodKids_denovo(allele1,allele1,allele1,allele2, i)*parentGLF[i][1];
  parentConditional[i][2] = likelihoodKids_denovo(allele1,allele1,allele2,allele2, i)*parentGLF[i][2];
  parentConditional[i][3] = likelihoodKids_denovo(allele1,allele2,allele1,allele1, i)*parentGLF[i][3];
  parentConditional[i][4] = likelihoodKids_denovo(allele1,allele2,allele1,allele2, i)*parentGLF[i][4];
  parentConditional[i][5] = likelihoodKids_denovo(allele1,allele2,allele2,allele2, i)*parentGLF[i][5];
  parentConditional[i][6] = likelihoodKids_denovo(allele2,allele2,allele1,allele1, i)*parentGLF[i][6];
  parentConditional[i][7] = likelihoodKids_denovo(allele2,allele2,allele1,allele2, i)*parentGLF[i][7];
  parentConditional[i][8] = likelihoodKids_denovo(allele2,allele2,allele2,allele2, i)*parentGLF[i][8];

  //Mariginal is joint likelihood of parents marginalizing all kids
  //This is to avoid redundant calculation
  parentMarginal[i][0] = parentConditional[i][0]*parentPrior[0];
  parentMarginal[i][1] = parentConditional[i][1]*parentPrior[1];
  parentMarginal[i][2] = parentConditional[i][2]*parentPrior[2];
  parentMarginal[i][3] = parentConditional[i][3]*parentPrior[3];
  parentMarginal[i][4] = parentConditional[i][4]*parentPrior[4];
  parentMarginal[i][5] = parentConditional[i][5]*parentPrior[5];
  parentMarginal[i][6] = parentConditional[i][6]*parentPrior[6];
  parentMarginal[i][7] = parentConditional[i][7]*parentPrior[7];
  parentMarginal[i][8] = parentConditional[i][8]*parentPrior[8];

}

void NucFamGenotypeLikelihood::CalcParentMarginal_leaveone(int i, double freq, int leave)
{
  double lkF11, lkF12, lkF22, lkM11, lkM12, lkM22; // lkC11, lkC12, lkC22;
  lkF11=lkF12=lkF22=lkM11=lkM12=lkM22=1.0;

  getGenoLikelihood(pedGLF->glf[i], allele1, allele2, &lkF11, &lkF12, &lkF22);
  getGenoLikelihood(pedGLF->glf[i]+1, allele1, allele2, &lkM11, &lkM12, &lkM22);

  parentGLF[i][0] = lkF11*lkM11;
  parentGLF[i][1] = lkF11*lkM12;
  parentGLF[i][2] = lkF11*lkM22;
  parentGLF[i][3] = lkF12*lkM11;
  parentGLF[i][4] = lkF12*lkM12;
  parentGLF[i][5] = lkF12*lkM22;
  parentGLF[i][6] = lkF22*lkM11;
  parentGLF[i][7] = lkF22*lkM12;
  parentGLF[i][8] = lkF22*lkM22;


  if(nFam>1)
    SetParentPrior(freq);
  else 
    SetParentPriorSingleTrio();

  parentConditional[i][0] = likelihoodKids_leaveone(allele1,allele1,allele1,allele1, i, leave)*parentGLF[i][0];
  parentConditional[i][1] = likelihoodKids_leaveone(allele1,allele1,allele1,allele2, i, leave)*parentGLF[i][1];
  parentConditional[i][2] = likelihoodKids_leaveone(allele1,allele1,allele2,allele2, i, leave)*parentGLF[i][2];
  parentConditional[i][3] = likelihoodKids_leaveone(allele1,allele2,allele1,allele1, i, leave)*parentGLF[i][3];
  parentConditional[i][4] = likelihoodKids_leaveone(allele1,allele2,allele1,allele2, i, leave)*parentGLF[i][4];
  parentConditional[i][5] = likelihoodKids_leaveone(allele1,allele2,allele2,allele2, i, leave)*parentGLF[i][5];
  parentConditional[i][6] = likelihoodKids_leaveone(allele2,allele2,allele1,allele1, i, leave)*parentGLF[i][6];
  parentConditional[i][7] = likelihoodKids_leaveone(allele2,allele2,allele1,allele2, i, leave)*parentGLF[i][7];
  parentConditional[i][8] = likelihoodKids_leaveone(allele2,allele2,allele2,allele2, i, leave)*parentGLF[i][8];

  //Mariginal is joint likelihood of parents marginalizing all kids
  //This is to avoid redundant calculation
  parentMarginal[i][0] = parentConditional[i][0]*parentPrior[0];
  parentMarginal[i][1] = parentConditional[i][1]*parentPrior[1];
  parentMarginal[i][2] = parentConditional[i][2]*parentPrior[2];
  parentMarginal[i][3] = parentConditional[i][3]*parentPrior[3];
  parentMarginal[i][4] = parentConditional[i][4]*parentPrior[4];
  parentMarginal[i][5] = parentConditional[i][5]*parentPrior[5];
  parentMarginal[i][6] = parentConditional[i][6]*parentPrior[6];
  parentMarginal[i][7] = parentConditional[i][7]*parentPrior[7];
  parentMarginal[i][8] = parentConditional[i][8]*parentPrior[8];

}

//likelihood of reads of all kids conditional on parental genotypes
double NucFamGenotypeLikelihood::likelihoodKids(int fa1, int fa2, int ma1, int ma2, int famIdx)
{
  double lkKids = 1.0;
  double lk=1.0;
  double lk11, lk12, lk22; 
  int famSize = pedGLF->ped->families[famIdx]->count;
  glfHandler ** glf = pedGLF->glf;
  for(int i=2; i<famSize; i++)
    {
      lk = likelihoodONEKid(&glf[famIdx][i], fa1, fa2, ma1, ma2);
      lkKids *= lk;
    }
  return(lkKids);
}

//likelihood of reads of ONE kid conditional on parental genotypes
double NucFamGenotypeLikelihood::likelihoodONEKid(glfHandler *glf, int fa1, int fa2, int ma1, int ma2)
{
  double lk=1.0;
  double lk11, lk12, lk22; 

  getGenoLikelihood(glf, allele1, allele2, &lk11, &lk12, &lk22);

      // father 1/1
      if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele1)
	{lk = lk11;      /* kidsCondLikelihood[famIdx][i][0] = lk;*/ }
      else if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele2) 
	{lk = 0.5 * (lk11+lk12);      /* kidsCondLikelihood[famIdx][i][1] = lk; */}
      else if(fa1==allele1 && fa2==allele1 && ma1==allele2 && ma2==allele2) 
	{lk = lk12;      /* kidsCondLikelihood[famIdx][i][2] = lk; */ }
      
      // father 1/2
      else if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele1)
	{lk = 0.5 * (lk11 + lk12);     /* kidsCondLikelihood[famIdx][i][3] = lk; */}
      else if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele2) 
	{lk = 0.25*lk11 + 0.5*lk12 + 0.25*lk22;       /*kidsCondLikelihood[famIdx][i][4] = lk;*/}
      else if(fa1==allele1 && fa2==allele2 && ma1==allele2 && ma2==allele2 )
	{lk = 0.5 * (lk12 + lk22);       /*kidsCondLikelihood[famIdx][i][5] = lk; */}
      
      // father 2/2
      else if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele1) 
	{lk = lk12;      /*kidsCondLikelihood[famIdx][i][6] = lk; */ }
      else if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele2)
	{lk = 0.5 * (lk12 + lk22);       /*kidsCondLikelihood[famIdx][i][7] = lk; */}
      else if(fa1==allele2 && fa2==allele2 && ma1==allele2 && ma2==allele2) 
	{lk = lk22;       /*kidsCondLikelihood[famIdx][i][8] = lk;*/ }
      
  return(lk);
}

double NucFamGenotypeLikelihood::likelihoodONEKid_denovo(glfHandler *glf, int fa1, int fa2, int ma1, int ma2)
{
    double lk = 1.0;
    const double *ptrlk = glf->GetLikelihoods(pedGLF->currentPos);

      // father 1/1
      if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele1)
	{  lk = CalcDenovoMutLk(ptrlk, allele1, allele1);      /* kidsCondLikelihood[famIdx][i][0] = lk;*/ }
      else if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele2) 
	{  lk = 0.5 * (CalcDenovoMutLk(ptrlk, allele1, allele1) + CalcDenovoMutLk(ptrlk, allele1, allele2) );      /* kidsCondLikelihood[famIdx][i][1] = lk; */}
      else if(fa1==allele1 && fa2==allele1 && ma1==allele2 && ma2==allele2) 
	{  lk = CalcDenovoMutLk(ptrlk, allele1, allele2);      /* kidsCondLikelihood[famIdx][i][2] = lk; */ }
      
      // father 1/2
      else if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele1)
	{  lk = 0.5 * (CalcDenovoMutLk(ptrlk, allele1, allele1) + CalcDenovoMutLk(ptrlk, allele1, allele2));     /* kidsCondLikelihood[famIdx][i][3] = lk; */}
      else if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele2) 
	{  lk = 0.25*CalcDenovoMutLk(ptrlk, allele1, allele1) + 0.5*CalcDenovoMutLk(ptrlk, allele1, allele2) + 0.25*CalcDenovoMutLk(ptrlk, allele2, allele2);       /*kidsCondLikelihood[famIdx][i][4] = lk;*/}
      else if(fa1==allele1 && fa2==allele2 && ma1==allele2 && ma2==allele2 )
	{  lk = 0.5 * (CalcDenovoMutLk(ptrlk, allele1, allele2) + CalcDenovoMutLk(ptrlk, allele2, allele2));       /*kidsCondLikelihood[famIdx][i][5] = lk; */}
      
      // father 2/2
      else if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele1) 
	{lk = CalcDenovoMutLk(ptrlk, allele1, allele2);      /*kidsCondLikelihood[famIdx][i][6] = lk; */ }
      else if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele2)
	{lk = 0.5 * (CalcDenovoMutLk(ptrlk, allele1, allele2) + CalcDenovoMutLk(ptrlk, allele2, allele2));       /*kidsCondLikelihood[famIdx][i][7] = lk; */}
      else if(fa1==allele2 && fa2==allele2 && ma1==allele2 && ma2==allele2) 
	{lk = CalcDenovoMutLk(ptrlk, allele2, allele2);       /*kidsCondLikelihood[famIdx][i][8] = lk;*/ }
	
	return(lk);
}

//likelihood of reads of all kids conditional on parental genotypes allowing for de novo mutations
double NucFamGenotypeLikelihood::likelihoodKids_denovo(int fa1, int fa2, int ma1, int ma2, int famIdx)
{
  double lkKids = 1.0;
  double lk=1.0;
  double lk11, lk12, lk22; 
  int famSize = pedGLF->ped->families[famIdx]->count;
  glfHandler ** glf = pedGLF->glf;
  for(int i=2; i<famSize; i++)
    {
      lk = likelihoodONEKid_denovo(&glf[famIdx][i], fa1, fa2, ma1, ma2);
      lkKids *= lk;
    }
  return(lkKids);
}

//likelihood of reads of kids except the ith one conditional on parental genotypes
double NucFamGenotypeLikelihood::likelihoodKids_leaveone(int fa1, int fa2, int ma1, int ma2, int famIdx, int leave)
{
  double lkKids = 1.0;
  double lk=1.0;
  double lk11, lk12, lk22; 
  int famSize = pedGLF->ped->families[famIdx]->count;
  glfHandler ** glf = pedGLF->glf;
  for(int i=2; i<famSize; i++)
    {
      if(i==leave) continue;
      lk = likelihoodONEKid(&glf[famIdx][i], fa1, fa2, ma1, ma2);
      lkKids *= lk;
    }
  return(lkKids);
}


//likelihood of the reads and a specific genotype of a kid conditional on parental genotypes
JointGenoLk NucFamGenotypeLikelihood::likelihoodKidGenotype(int fa1, int fa2, int ma1, int ma2, int famIdx, int kidIndex)
{
  double lkKidG11 = 1.0;
  double lkKidG12 = 1.0;
  double lkKidG22 = 1.0;
  double lk=0.0;
  double lk11, lk12, lk22;
  double lkg11, lkg12, lkg22;
  glfHandler **glf = pedGLF->glf;
  int famSize = pedGLF->ped->families[famIdx]->count;

  for(int i=2; i<famSize; i++)
    {
      getGenoLikelihood(&glf[famIdx][i], allele1, allele2, &lk11, &lk12, &lk22);
      
      // father 1/1
      if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele1)
	{lk = lk11; lkg11 = lk11; lkg12=lkg22=0; }
      if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele2) 
	{lk = 0.5 * (lk11+lk12); lkg11=lk11*0.5; lkg12=lk12*0.5; lkg22=0;}
      if(fa1==allele1 && fa2==allele1 && ma1==allele2 && ma2==allele2) 
	{lk = lk12; lkg11=0; lkg12=lk12; lkg22=0; }
      
      // father 1/2
      if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele1)
	{lk = 0.5 * (lk11 + lk12); lkg11=lk11*0.5; lkg12=lk12*0.5; lkg22=0; }
      if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele2) 
	{lk = 0.25*lk11 + 0.5*lk12 + 0.25*lk22; lkg11=lk11*0.25; lkg12=lk12*0.5; lkg22=lk22*0.25;}
      if(fa1==allele1 && fa2==allele2 && ma1==allele2 && ma2==allele2)
	{lk = 0.5 * (lk12 + lk22); lkg11=0; lkg12=lk12*0.5; lkg22=lk22*0.5;}
      
      // father 2/2
      if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele1) 
	{lk = lk12; lkg11=0; lkg12=lk12; lkg22=0; }
      if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele2)
	{lk = 0.5 * (lk12 + lk22); lkg11=0; lkg12=lk12*0.5; lkg22=lk22*0.5;}
      if(fa1==allele2 && fa2==allele2 && ma1==allele2 && ma2==allele2) 
	{lk = lk22; lkg11=0; lkg12=0; lkg22=lk22;}
        
      if(i!=kidIndex)
	{
	  lkKidG11 *= lk;
	  lkKidG12 *= lk;
	  lkKidG22 *= lk;
	}
      else
	{
	  lkKidG11*=lkg11;
	  lkKidG12*=lkg12;
	  lkKidG22*=lkg22;
	}
    }
  JointGenoLk JGLK;
  JGLK.g11 = lkKidG11;
  JGLK.g12 = lkKidG12;
  JGLK.g22 = lkKidG22;
  
  return(JGLK);
}

//likelihood of the reads and a specific genotype of a kid conditional on parental genotypes, allowing for de novo mutations
JointGenoLk_denovo NucFamGenotypeLikelihood::likelihoodKidGenotype_denovo(int fa1, int fa2, int ma1, int ma2, int famIdx, int kidIndex)
{
  double lk=1.0;
  glfHandler **glf = pedGLF->glf;
  int famSize = pedGLF->ped->families[famIdx]->count;
  double lkKidGeno[10];
  
  for(int i=0; i<10; i++) lkKidGeno[i] = 1.0;
  
  for(int i=2; i<famSize; i++)
    {
      if(i!=kidIndex)
	{
	  lk = likelihoodONEKid_denovo(&glf[famIdx][i], fa1, fa2, ma1, ma2);
	  for(int k=0; k<10; k++)
	    lkKidGeno[k] *= lk;
	}
      else
	{
	  double lkKidGenoJoint[10];
	  GetJointGenoLk_denovo(&glf[famIdx][i], fa1, fa2, ma1, ma2, lkKidGenoJoint);
	  for(int k=0; k<10; k++)
	    lkKidGeno[k] *= lkKidGenoJoint[k];
	}
    }
  
  JointGenoLk_denovo JGLK;
  
  for(int i=0; i<10; i++)
    JGLK.geno[i] = lkKidGeno[i];
  
  return(JGLK);
}

void NucFamGenotypeLikelihood::GetJointGenoLk_denovo(glfHandler *glf, int fa1, int fa2, int ma1, int ma2, double *lkKidGenoJoint)
{
  const double *ptrlk = glf->GetLikelihoods(pedGLF->currentPos);
  int idx1, idx2, idx3;
  
  // father 1/1
  if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele1)
    {
      idx1 = glfHandler::GenotypeIndex(allele1, allele1);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = gM.genoMutMatrix[idx1][i]*ptrlk[i];
    }
  if(fa1==allele1 && fa2==allele1 && ma1==allele1 && ma2==allele2) 
    {
      idx1 = glfHandler::GenotypeIndex(allele1, allele1);
      idx2 = glfHandler::GenotypeIndex(allele1, allele2);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = (0.5*gM.genoMutMatrix[idx1][i] + 0.5*gM.genoMutMatrix[idx2][i]) * ptrlk[i];
    }
  if(fa1==allele1 && fa2==allele1 && ma1==allele2 && ma2==allele2) 
    {
      idx2 = glfHandler::GenotypeIndex(allele1, allele2);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = gM.genoMutMatrix[idx2][i] * ptrlk[i];
    }
  
  // father 1/2
  if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele1)
    {
      idx1 = glfHandler::GenotypeIndex(allele1, allele1);
      idx2 = glfHandler::GenotypeIndex(allele1, allele2);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = (0.5*gM.genoMutMatrix[idx1][i] + 0.5*gM.genoMutMatrix[idx2][i]) * ptrlk[i];
    }
  if(fa1==allele1 && fa2==allele2 && ma1==allele1 && ma2==allele2) 
    {
      idx1 = glfHandler::GenotypeIndex(allele1, allele1);
      idx2 = glfHandler::GenotypeIndex(allele1, allele2);
      idx3 = glfHandler::GenotypeIndex(allele2, allele2);	 
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = (0.25*gM.genoMutMatrix[idx1][i] + 0.5*gM.genoMutMatrix[idx2][i] + 0.25*gM.genoMutMatrix[idx3][i]) * ptrlk[i];
    }
  if(fa1==allele1 && fa2==allele2 && ma1==allele2 && ma2==allele2)
    {
      idx2 = glfHandler::GenotypeIndex(allele1, allele2);
      idx3 = glfHandler::GenotypeIndex(allele2, allele2);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = (0.5*gM.genoMutMatrix[idx2][i] + 0.5*gM.genoMutMatrix[idx3][i]) * ptrlk[i];
    }
  
  // father 2/2
  if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele1) 
    {
      idx2 = glfHandler::GenotypeIndex(allele1, allele2);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = gM.genoMutMatrix[idx2][i] * ptrlk[i];
    }
  if(fa1==allele2 && fa2==allele2 && ma1==allele1 && ma2==allele2)
    {
      idx2 = glfHandler::GenotypeIndex(allele1, allele2);
      idx3 = glfHandler::GenotypeIndex(allele2, allele2);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = (0.5*gM.genoMutMatrix[idx2][i] + 0.5*gM.genoMutMatrix[idx3][i]) * ptrlk[i];	
    }
  if(fa1==allele2 && fa2==allele2 && ma1==allele2 && ma2==allele2) 
    {
      idx3 = glfHandler::GenotypeIndex(allele2, allele2);
      for(int i=0; i<10; i++)
	lkKidGenoJoint[i] = gM.genoMutMatrix[idx3][i] * ptrlk[i];		
    }
  
}

double NucFamGenotypeLikelihood::CalcDenovoMutLk(const double *ptrlk, int a1, int a2)
{
  double lk = 0.0;
  int idx = glfHandler::GenotypeIndex(a1, a2);
  
  for(int i=0; i<10; i++)
    lk += gM.genoMutMatrix[idx][i] * ptrlk[i];
  
  return(lk);
}

int NucFamGenotypeLikelihood::GetBestGenoIdx(double p11, double p12, double p22)
{
  int bestIdx = 0;
  double best = p11;
  if(p12>best) {best = p12; bestIdx = 1;}
  if(p22>best) {best = p22; bestIdx = 2; }
  return(bestIdx);
}

String NucFamGenotypeLikelihood::GetBestGenoLabel(int best){
  String genotypeLabel[10] = {"A/A", "A/C", "A/G", "A/T", "C/C", "C/G", "C/T", "G/G", "G/T", "T/T"};
  int labelIdx = 0;
  if(best==0)
    labelIdx = glfHandler::GenotypeIndex(allele1, allele1);
  if(best==1)
    labelIdx = glfHandler::GenotypeIndex(allele1, allele2);
  if(best==2)
    labelIdx = glfHandler::GenotypeIndex(allele2, allele2);
  return(genotypeLabel[labelIdx]);
}

String NucFamGenotypeLikelihood::GetBestGenoLabel_denovo(int best){
  String genotypeLabel[10] = {"A/A", "A/C", "A/G", "A/T", "C/C", "C/G", "C/T", "G/G", "G/T", "T/T"};
  return(genotypeLabel[best]);
}

String NucFamGenotypeLikelihood::GetBestGenoLabel_vcfv4(int best){
  String genotypeLabel[5] = {"0/0", "0/1", "1/1", "1/2", "2/2"};
  int labelIdx;
  if(best==0)
    if(pedGLF->refBase==allele1) labelIdx = 0; else labelIdx = 2;
  if(best==1)
    if(pedGLF->refBase==allele1) labelIdx = 1; else labelIdx = 3;
  if(best==2)
    if(pedGLF->refBase==allele1) labelIdx = 2 ; else labelIdx = 4; 

  return(genotypeLabel[labelIdx]);
}

double NucFamGenotypeLikelihood::CalcDosage(int i, int j){
  double dose = postProb[i][j][1] + postProb[i][j][2]*2;
  return(dose);
}

//get 3 conditional likelihood of reads given genotypes
void NucFamGenotypeLikelihood::getLogGenoLikelihood(glfHandler *glf,int a1, int a2, unsigned char *lk11, unsigned char *lk12, unsigned char *lk22)
{  
  if(glf->handle==NULL) 
  {
    *lk11 = 0.0;
    *lk12 = 0.0; 
    *lk22 = 0.0;
    return;
  }

  const unsigned char *ptrlk = glf->GetLogLikelihoods(pedGLF->currentPos);
 
  *lk11 = ptrlk[geno11];
  *lk12 = ptrlk[geno12];
  *lk22 = ptrlk[geno22];  
}

void NucFamGenotypeLikelihood::getGenoLikelihood(glfHandler *glf,int a1, int a2, double *lk11, double *lk12, double *lk22)
{
  if(glf->handle==NULL) 
    { 
      *lk11 = 1.0;
      *lk12 = 1.0; 
      *lk22 = 1.0;
      return;
    }
  
  const double *ptrlk = glf->GetLikelihoods(pedGLF->currentPos);
  
  *lk11 = ptrlk[geno11];
  *lk12 = ptrlk[geno12];
  *lk22 = ptrlk[geno22];
}

int NucFamGenotypeLikelihood::CalcMaxLogLkIdx(std::vector<double> &loglk, int n)
{
  int idx = 0;
  double max = loglk[0];
  
  for(int i=0; i<n; i++){
    if(max<loglk[i]) {
      max = loglk[i];
      idx = i;
    }
  }
  return(idx);
}

double NucFamGenotypeLikelihood::CalcSum(std::vector<double> &ratio)
{
  double sum = 0.0;
  for(uint i=0; i<ratio.size(); i++)
    sum += ratio[i];
  return(sum);
}

int NucFamGenotypeLikelihood::CalcVarPosterior(int n)
{
  int allele1, allele2;
  int ts       = (Poly::ts(pedGLF->refBase));   //transition
  int tvs1     = (Poly::tvs1(pedGLF->refBase)); //transvertions1
  int tvs2     = (Poly::tvs2(pedGLF->refBase)); //transvertion2
  
  
  int maxidx = CalcMaxLogLkIdx(varllk, n);
  std::vector<double> ratio; 
  ratio.resize(n);
  
  for(uint i=0; i<ratio.size(); i++){
    ratio[i] = pow10(varllk[i]-varllk[maxidx]);
  }
  
  double sumRatio = CalcSum(ratio);
  varPostProb = 1/sumRatio;
  
  if(varPostProb<par->posterior) return(-1);
  
  if(maxidx==0){
    allele1 = allele2 = pedGLF->refBase; 
  }
  else if(maxidx==1){
    allele1 = pedGLF->refBase; allele2 = ts; 
  }
  else if(maxidx==2){ 
    allele1 = pedGLF->refBase; allele2 = tvs1; 
  }
  else if(maxidx==3){ 
    allele1 = pedGLF->refBase; allele2 = tvs2; 
  }
  else if(maxidx==4){ 
    allele1 = ts; allele2 = tvs1;  
  }
  else if(maxidx==5) { 
    allele1 = ts; allele2 = tvs2;
  }
  else if(maxidx==6){
    allele1 = tvs1; allele2 = tvs2; 
  }
  
  SetAlleles(allele1, allele2);
  SetMleFreq(varfreq[maxidx]);
  SetParentPrior(varfreq[maxidx]);
  CalcPolyQual(varPostProb);

  return(maxidx);
}

void NucFamGenotypeLikelihood::CalcPolyQual(double postProb)
{
  double qual;
  if(postProb>0.9999999999) qual = 100;
  else qual = -10*log10(1-postProb); 
  SetPolyQual(qual);
}

void NucFamGenotypeLikelihood::OutputVCF(FILE * fh)
{
  if(!vcf) {
    time_t t; time(&t);
    fprintf(fh, "##fileformat=VCFv4.0\n");
    fprintf(fh, "##fileDate=%s", ctime(&t));
    fprintf(fh, "##minMapQuality=%f\n", par->minMapQuality);
    fprintf(fh, "##minTotalDepth=%d\n", par->minTotalDepth);
    fprintf(fh, "##maxTodalDepth=%d\n", par->maxTotalDepth);
    fprintf(fh, "##posterior=%.3f\n",   par->posterior);
    fprintf(fh, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
    fprintf(fh, "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Percentage of Samples With Data\">\n");
    fprintf(fh, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n");
    if(nFam>1 || ped->families[0]->isNuclear()==false) 
    fprintf(fh, "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Reference Allele Frequency\">\n");
    fprintf(fh, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(fh, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
    fprintf(fh, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    fprintf(fh, "##FORMAT=<ID=DS,Number=1,Type=Float, Description=\"Dosage: Defined As the Expected Alternative Allele Count\">\n");
    if(!par->gl_off) fprintf(fh, "##FORMAT=<ID=GL,Number=3,Type=Unsigned Char, Description=\"Genotype Likelihoods\">\n");
    fprintf(fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for(int i=0; i<nFam; i++)
      for(int j=0; j<pedGLF->ped->families[i]->count; j++)
	{
	  int idx = pedGLF->ped->families[i]->path[j];
	  String famid = pedGLF->ped->families[i]->famid;
	  String pid =  pedGLF->ped->persons[idx]->pid;
	  fprintf(fh, "\t%s", (const char*)pid);
	}
    fprintf(fh, "\n");
    vcf = true;
  }
 
  String FILTER;
  FILTER.printf(".");

  String INFO;
  if(nFam==1 && ped->families[0]->isNuclear())
   INFO.printf("NS=%d;PS=%.1f;DP=%d", numSampWithData, percSampWithData*100, totalDepth);
  else
  INFO.printf("NS=%d;PS=%.1f;DP=%d;AF=%.4f", numSampWithData, percSampWithData*100, totalDepth, GetMinimizer());

  String FORMAT = "GT:GQ:DP:DS"; if(!par->gl_off) FORMAT+=":GL";
  char bases[5]  = {'0', 'A', 'C', 'G', 'T'};
  String alt;
  if(pedGLF->refBase==allele1) alt = bases[allele2]; else alt = alt+bases[allele1]+","+bases[allele2];

  fprintf(fh, "%s\t%d\t%s\t%c\t%s\t%d\t%s\t%s\t%s", pedGLF->GetNonNULLglf()->label.c_str(), pedGLF->currentPos+1, ".", bases[pedGLF->refBase], alt.c_str(), int(polyQual+0.5), FILTER.c_str(), INFO.c_str(), FORMAT.c_str());
 
  int GTQual = 0.0;
  int bestIdx = 0;
  unsigned char llk11, llk12, llk22;
  for(int i=0; i<nFam; i++)
    for(int j=0; j<pedGLF->ped->families[i]->count; j++)
      {
	bestIdx = bestGenoIdx[i][j];
	if(postProb[i][j][bestIdx]>0.9999999999) GTQual = 100;
	else GTQual = int(-10.*log10(1.-postProb[i][j][bestIdx])+0.5);
	getLogGenoLikelihood(&pedGLF->glf[i][j], allele1, allele2, &llk11, &llk12, &llk22);
	fprintf(fh, "\t%s:", bestGenoLabel[i][j].c_str());
	fprintf(fh, "%d:", GTQual);
	fprintf(fh, "%d:", pedGLF->glf[i][j].handle==NULL? 0 : pedGLF->glf[i][j].GetDepth(pedGLF->currentPos));
	fprintf(fh, "%.2f", dosage[i][j]);
	if(!par->gl_off)fprintf(fh, ":%u,%u,%u", llk11, llk12, llk22);

      }
  fprintf(fh, "\n");
  fflush(fh);
}

void NucFamGenotypeLikelihood::OutputVCF_denovo(FILE * fh)
{
  if(!vcf) {
    time_t t; time(&t);
    fprintf(fh, "##fileformat=VCFv4.0\n");
    fprintf(fh, "##fileDate=%s", ctime(&t));
    fprintf(fh, "##minMapQuality=%f\n", par->minMapQuality);
    fprintf(fh, "##minTotalDepth=%d\n", par->minTotalDepth);
    fprintf(fh, "##maxTodalDepth=%d\n", par->maxTotalDepth);
    fprintf(fh, "##posterior=%.3f\n",   par->posterior);
    fprintf(fh, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
    fprintf(fh, "##INFO=<ID=PS,Number=1,Type=Integer,Description=\"Percentage of Samples With Data\">\n");
    fprintf(fh, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">\n");
    if(nFam>1 || ped->families[0]->isNuclear()==false) 
    fprintf(fh, "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Reference Allele Frequency\">\n");
    fprintf(fh, "##INFO=<ID=DQ,Number=1,Type=Float, Description=\"De Novo Mutation Quality\">\n");
    fprintf(fh, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(fh, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
    fprintf(fh, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    if(!par->gl_off) fprintf(fh, "##FORMAT=<ID=GL,Number=10,Type=Unsigned Char, Description=\"Genotype Likelihoods\">\n");
    fprintf(fh, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    for(int i=0; i<nFam; i++)
      for(int j=0; j<pedGLF->ped->families[i]->count; j++)
	{
	  int idx = pedGLF->ped->families[i]->path[j];
	  String famid = pedGLF->ped->families[i]->famid;
	  String pid =  pedGLF->ped->persons[idx]->pid;
	  fprintf(fh, "\t%s", (const char*)pid);
	}
    fprintf(fh, "\n");
    vcf = true;
  }

  if(denovoLR < par->denovoLR) return;

  if(denovo_mono) allele2=allele1;
  denovo_mono = false;

  String FILTER;
  FILTER.printf(".");

  String INFO;

  if(nFam==1 && ped->families[0]->isNuclear())
   INFO.printf("NS=%d;PS=%.1f;DP=%d;DQ=%.3f", numSampWithData, percSampWithData*100, totalDepth, denovoLR);
  else
   INFO.printf("NS=%d;PS=%.1f;DP=%d;AF=%.4f;DQ=%.3f", numSampWithData, percSampWithData*100, totalDepth, GetMinimizer(), denovoLR);

  String FORMAT = "GT:GQ:DP"; if(par->gl_off==false) FORMAT+=":GL";
  char bases[5]  = {'0', 'A', 'C', 'G', 'T'};
  String alt;
  if(pedGLF->refBase==allele1) alt = bases[allele2]; else alt = alt+bases[allele1]+","+bases[allele2];

  fprintf(fh, "%s\t%d\t%s\t%c\t%s\t%d\t%s\t%s\t%s", pedGLF->GetNonNULLglf()->label.c_str(), pedGLF->currentPos+1, ".", bases[pedGLF->refBase], alt.c_str(), int(polyQual+0.5), FILTER.c_str(), INFO.c_str(), FORMAT.c_str());

  int GTQual = 0.0;
  int bestIdx = 0;
  const unsigned char *ptrlk;

  for(int i=0; i<nFam; i++)
    for(int j=0; j<pedGLF->ped->families[i]->count; j++)
      {
	ptrlk = pedGLF->glf[i][j].GetLogLikelihoods(pedGLF->currentPos);
	bestIdx = bestGenoIdx[i][j];
	if(postProb[i][j][bestIdx]>0.9999999999) GTQual = 100;
	else GTQual = int(-10.*log10(1.-postProb[i][j][bestIdx])+0.5);
	fprintf(fh, "\t%s:", bestGenoLabel[i][j].c_str());
	fprintf(fh, "%d:", GTQual);
	fprintf(fh, "%d", pedGLF->glf[i][j].handle==NULL? 0 : pedGLF->glf[i][j].GetDepth(pedGLF->currentPos));

	if(par->gl_off==false)
 	{
	//Output 10 GL values
	fprintf(fh, ":");
	for(int g=0; g<9; g++)
	  fprintf(fh, "%d,", ptrlk[g]);
	fprintf(fh, "%d",ptrlk[9]);
	}
      }
  fprintf(fh, "\n");
  fflush(fh);
}

