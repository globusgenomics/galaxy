#ifndef __MUTATIONMODEL_H__
#define __MUTATIONMODEL_H__
#include "StringArray.h"
#include <math.h>

class AlleleMutationModel
{
 private:
 double mu;
 double tstv;
 
 public:
 double alleleMutMatrix[4][4];
 
 public:
 AlleleMutationModel();
 ~AlleleMutationModel();
 void SetAlleleMutMatrix();
 void SetMutationRate(double);
 void SetTsTvRatio(double);
 void PrintMatrix();
};

class GenotypeMutationModel
{
 public:
 GenotypeMutationModel();
 ~GenotypeMutationModel();
 
 double genoMutMatrix[10][10];
 void SetGenoMutMatrix(AlleleMutationModel&);
 void PrintMatrix();
};

#endif
