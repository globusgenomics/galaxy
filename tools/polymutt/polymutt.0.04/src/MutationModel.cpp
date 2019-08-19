#include "MutationModel.h"

AlleleMutationModel::AlleleMutationModel()
{
 mu=0;
 tstv = 0.5;
}

AlleleMutationModel::~AlleleMutationModel(){}

void AlleleMutationModel::SetMutationRate(double r){mu=r;}

void AlleleMutationModel::SetTsTvRatio(double r) {tstv = r; }

void AlleleMutationModel::SetAlleleMutMatrix()
{
  for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
      if(i==j)
	alleleMutMatrix[i][j] = 1-mu;
      else
	alleleMutMatrix[i][j] = (1-mu)/3;
  
  if(tstv!=0.0)
    {
      alleleMutMatrix[0][2] = alleleMutMatrix[2][0] = alleleMutMatrix[1][3] = alleleMutMatrix[3][1] = mu/3*(3-3/(1+tstv));
      alleleMutMatrix[0][1] = alleleMutMatrix[0][3] = alleleMutMatrix[1][0] = alleleMutMatrix[1][2] =
	alleleMutMatrix[2][1] = alleleMutMatrix[2][3] = alleleMutMatrix[3][0] = alleleMutMatrix[3][2] = mu/3*(0.5/(1+tstv)*3);
    }
}

void AlleleMutationModel::PrintMatrix()
{
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<4; j++)
	printf("%f ", alleleMutMatrix[i][j]);
      printf("\n");
    }
}


GenotypeMutationModel::GenotypeMutationModel(){}
GenotypeMutationModel::~GenotypeMutationModel(){}

void GenotypeMutationModel::SetGenoMutMatrix(AlleleMutationModel &aM)
{
  double mutRate16[16][16];
  double mutSum = 0.0;
  int fromIdx = -1;
  int toIdx = -1;
  
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<4; j++)
	{
	  fromIdx++;
	  toIdx = -1;
	  mutSum = 0;
	  
	  for(int ii=0; ii<4; ii++)
	    {
	      for(int jj=0; jj<4; jj++)
		{ 
		  toIdx++; 
		  mutRate16[fromIdx][toIdx] = (aM.alleleMutMatrix[i][ii]*aM.alleleMutMatrix[j][jj]); 
		}
	    }
	}
    }
  
  
  // sum ordered heterozygotes transition probabilities
  int hetOrdered1[6] = {2,3,4,7,8,12};
  int hetOrdered2[6] = {5,9,13,10,14,15};
  
  for(int i=0; i<6; i++)
    { 
      for(int j=0; j<16; j++)
	{ 
	  mutRate16[j][hetOrdered1[i]-1] += mutRate16[j][hetOrdered2[i]-1];
	}
    }
  
  int unOrderedIdx[10] = {1,2,3,4,6,7,8,11,12,16};
  
  for(int i=0; i<10; i++)
    for(int j=0; j<10; j++)
      genoMutMatrix[i][j] = mutRate16[unOrderedIdx[i]-1][unOrderedIdx[j]-1];
}

void GenotypeMutationModel::PrintMatrix()
{
  for(int i=0; i<10; i++)
    {
      for(int j=0; j<10; j++)
	printf("%f ", genoMutMatrix[i][j]);
      printf("\n");
    }
}
