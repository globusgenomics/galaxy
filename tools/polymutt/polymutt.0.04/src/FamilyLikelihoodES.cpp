#include "FamilyLikelihoodES.h"
#include <map>
#include <vector>
#include "Error.h"

using namespace std;

ES_Peeling::ES_Peeling()
{
  ped = NULL;
  family = NULL;
  famIdx = -1;
  parents = NULL;
  offspring = NULL;
  spouses = NULL;
}

ES_Peeling::~ES_Peeling()
{
  if(parents==NULL) delete[] parents;
  if(offspring!=NULL) delete[] offspring;
  if(spouses!=NULL) delete[] spouses;
}

void ES_Peeling::SetPedigree(Pedigree *pedptr){ped = pedptr;}

void ES_Peeling::SetFamily(Family *famptr)
{
  if(family)
    {
      delete[] parents;
      delete[] offspring;
      delete[] spouses;
    }
  family = famptr;
  ped = &family->ped;
  famSize = family->count;
  parents = new IntArray[famSize];
  offspring = new IntArray[famSize];
  spouses = new IntArray[famSize];
}

void ES_Peeling::SetFamily(int fam){famIdx = fam;}


void ES_Peeling::SetupConnections()
{
  Person *p;
  int fa_traverse, mo_traverse;
  
  std::map<std::pair<int, int>, int> couples;
  std::pair<int, int> couple; 
  
  for(int i=0; i<famSize; i++)
    {
      p = ped->persons[family->path[i]];
      if(p->isFounder()) continue;
      
      fa_traverse = ped->persons[family->path[i]]->father->traverse;
      mo_traverse = ped->persons[family->path[i]]->mother->traverse;
      
      parents[i].Push(fa_traverse);
      parents[i].Push(mo_traverse);
      
      offspring[fa_traverse].Push(i);
      offspring[mo_traverse].Push(i);
      
      couple.first = fa_traverse;
      couple.second = mo_traverse;  
      
      if(couples[couple]==0)
	{
	  spouses[fa_traverse].Push(mo_traverse);
	  spouses[mo_traverse].Push(fa_traverse);
	}
      couples[couple]++;
    }
}

void ES_Peeling::BuildInitialPeelable()
{
  std::pair<int, int> couple;
  std::map<int, int> roof_visited;
  for(int i=0; i<famSize; i++)
    { 
      if(isLeaf(i))
	{
	  leaf.Push(i);
	  continue;
	}
      
      if(isRoof(i))
	{
	  if(roof_visited[i]>0 || roof_visited[spouses[i][0]]>0) continue;
	  if(ped->persons[family->path[i]]->sex==1)
	    {
	      couple.first = i;
	      couple.second = spouses[i][0];
	    }
	  else {
	    couple.second = i;
	    couple.first = spouses[i][0];
	  }
	  roof.push_back(couple); 
	  roof_visited[i]++; roof_visited[spouses[i][0]]++;
	  continue;
	}
      
      if(isPeripheral(i))
	{
	  peripheral.Push(i);
	  continue;
	} 
    }
}

void ES_Peeling::PrintPeelable()
{
  printf("Leaf:");
  for(int i=0; i<leaf.Length(); i++)
    printf(" %d", leaf[i]);
  printf("\n");
  
  printf("Peripheral:");
  for(int i=0; i<peripheral.Length(); i++)
    printf(" %d", peripheral[i]);
  printf("\n");
  
  printf("Roof:");
  for(int i=0; i<roof.size(); i++)
    printf(" %d-%d", roof[i].first, roof[i].second);
  printf("\n");
}

void ES_Peeling::BuildPeelingOrder()
{
  std::pair<int, int> peelFrom;
  std::pair<int, int> peelTo;
  int aLeaf, aPeripheral;
  std::pair<int, int> aRoof;
  std::vector<int> peeled;
  bool done = false;
  
  for(;;)
    {
      if(leaf.Length()==0 && roof.size()==0 && peripheral.Length()==0) break;
      
      if(done) break;
      
      while(leaf.Length()>0)
	{
	  aLeaf = leaf[0];
	  leaf.Delete(0);
	  
	  peeled.push_back(aLeaf);
	  
	  peelFrom.first = aLeaf;
	  peelFrom.second = -1;
	  from.push_back(peelFrom);
	  peelingType.push_back(1);
	  
	  peelTo.first = parents[aLeaf][0]; 
	  peelTo.second = parents[aLeaf][1];
	  to.push_back(peelTo);
	  
	  int idx = offspring[peelTo.first].Find(aLeaf);
	  if(idx<0) error("Peeling error for person %s in family %s! Check pedigree structure!!\n", ped->persons[family->path[aLeaf]]->pid.c_str(), family->famid.c_str());
	  offspring[peelTo.first].Delete(idx);
	  
	  idx = offspring[peelTo.second].Find(aLeaf);
	  if(idx<0) error("Peeling leaf error: %d!\n", aLeaf);
	  offspring[peelTo.second].Delete(idx);
	  parents[aLeaf].Delete(0);
	  parents[aLeaf].Delete(0);
	  
	  if(isPeripheral(peelTo.first)) 
	    { 
	      peripheral.Push(peelTo.first); 
	    }
	  if(isPeripheral(peelTo.second))
	    {
	      peripheral.Push(peelTo.second);
	    }
	  
	  int pos = find_element(roof, peelTo);
	  if(pos>0)
	    remove_element(roof, pos);
	  
	  if(peeled.size()==famSize-1) done = true;;
	}
      
      if(done) break;
      
      while(peripheral.Length()>0)
	{
	  aPeripheral = peripheral[0];
	  peripheral.Delete(0);
	  
	  peeled.push_back(aPeripheral);
	  
	  peelFrom.first = aPeripheral;
	  peelFrom.second = -1;
	  from.push_back(peelFrom);
	  peelingType.push_back(2);
	  
	  peelTo.first = spouses[aPeripheral][0];
	  peelTo.second = -1;
	  to.push_back(peelTo);
	  
	  if(spouses[aPeripheral].Length()>1)
	    error("Peripheral parent can not have more than one spouses!\n");
	  
	  int idx = spouses[spouses[aPeripheral][0]].Find(aPeripheral);
	  if(idx<0) error("No spouse can be found for person with PID of %s %d!\n", (const char *) ped->persons[family->path[aPeripheral]]->pid);
	  
	  spouses[spouses[aPeripheral][0]].Delete(idx);
	  spouses[aPeripheral].Delete(0);
	  
	  if(isFinal(peelTo.first))
	    {
	      if(peeled.size()!=famSize-1)
		error("Inconsistency observed!\n");
	      done = true;
	      break;
	    }
	  
	  if(isLeaf(spouses[aPeripheral][0]))
	    leaf.Push(spouses[aPeripheral][0]);
	  else if(isPeripheral(spouses[aPeripheral][0]))
	    peripheral.Push(spouses[aPeripheral][0]);
	  else if(isRoof(spouses[aPeripheral][0]))
	    UpdateRoof(roof, spouses[aPeripheral][0]);     
	}
      
      if(done) break;
      
      if(leaf.Length()>0) continue;
      if(peripheral.Length()>0) continue;
      
      while(roof.size()>0)
	{
	  aRoof = roof[0];
	  roof.erase(roof.begin()+0);
	  if(offspring[aRoof.first].Length()!=1 || offspring[aRoof.second].Length()!=1)
	    error("Roof can only have one offspring for peeling!\n");
	  
	  peeled.push_back(aRoof.first);
	  peeled.push_back(aRoof.second);
	  
	  from.push_back(aRoof);
	  peelingType.push_back(3);
	  
	  peelTo.first = offspring[aRoof.first][0];
	  peelTo.second = -1;
	  to.push_back(peelTo);
	  
	  parents[peelTo.first].Delete(0);
	  parents[peelTo.first].Delete(0);
	  offspring[aRoof.first].Delete(0);
	  offspring[aRoof.second].Delete(0);
	  
	  if(isPeripheral(peelTo.first))
	    peripheral.Push(peelTo.first);
	  else if(isRoof(peelTo.first))
	    UpdateRoof(roof, peelTo.first);
	  else if(isFinal(peelTo.first))
	    {
	      done = true;
	      break;
	    }  
	}
      if(done) break;
    }

   if(peeled.size()< famSize-1)
        error("Are there inbreeding loops in the pedigree? It cannot handel inbreeding yet!\n");

}

void ES_Peeling::BuildInitialPeelable2()
{
  std::pair<int, int> couple;
  std::map<int, int> roof_visited;
  
  for(int i=famSize-1; i>=0; i--)
    {       
      if(isRoof(i))
	{
	  if(roof_visited[i]>0 || roof_visited[spouses[i][0]]>0) continue;
	  if(ped->persons[family->path[i]]->sex==1)
	    {
	      couple.first = i;
	      couple.second = spouses[i][0];
	    }
	  else {
	    couple.second = i;
	    couple.first = spouses[i][0];
	  }
	  roof.push_back(couple); 
	  roof_visited[i]++; roof_visited[spouses[i][0]]++;
	  continue;
	}
      
      if(isPeripheral(i))
	{
	  peripheral.Push(i);
	  continue;
	} 
      
      if(isLeaf(i))
	{
	  leaf.Push(i);
	  continue;
	}
    }
}



void ES_Peeling::BuildPeelingOrder2()
{
  std::pair<int, int> peelFrom;
  std::pair<int, int> peelTo;
  int aLeaf, aPeripheral;
  std::pair<int, int> aRoof;
  std::vector<int> peeled;
  bool done = false;
  PrintPeelable(); 
  for(;;)
    {
      if(leaf.Length()==0 && roof.size()==0 && peripheral.Length()==0) break;
      
      if(roof.size()>0)
	{
	  aRoof = roof[0];
	  roof.erase(roof.begin()+0);
	  if(offspring[aRoof.first].Length()!=1 || offspring[aRoof.second].Length()!=1)
	    error("Roof can only have one offspring for peeling: %s %s!\n", GetPID(aRoof.first).c_str(), GetPID(aRoof.second).c_str());
	  
	  peeled.push_back(aRoof.first);
	  peeled.push_back(aRoof.second);
          
	  from.push_back(aRoof);
	  peelingType.push_back(3);
          
	  peelTo.first = offspring[aRoof.first][0];
	  peelTo.second = -1;
	  to.push_back(peelTo);
          
	  parents[peelTo.first].Delete(0);
	  parents[peelTo.first].Delete(0);
	  offspring[aRoof.first].Delete(0);
	  offspring[aRoof.second].Delete(0);
          
	  if(isPeripheral(peelTo.first))
	    peripheral.Push(peelTo.first);
	  else if(isRoof(peelTo.first))
	    UpdateRoof(roof, peelTo.first);
	  else if(isFinal(peelTo.first))
	    break;
	  continue;
	}
      
      if(peripheral.Length()>0)
	{
	  aPeripheral = peripheral[0];
	  peripheral.Delete(0);
          
	  peeled.push_back(aPeripheral);
          
	  peelFrom.first = aPeripheral;
	  peelFrom.second = -1;
	  from.push_back(peelFrom);
	  peelingType.push_back(2);
          
	  peelTo.first = spouses[aPeripheral][0];
	  peelTo.second = -1;
	  to.push_back(peelTo);
          
	  if(spouses[aPeripheral].Length()>1)
	    error("Peripheral parent can not have more than one spouses!\n");
	  
	  int idx = spouses[spouses[aPeripheral][0]].Find(aPeripheral);
	  if(idx<0) error("No spouse can be found for person with PID of %s %d!\n", (const char *) ped->persons[family->path[aPeripheral]]->pid);
          
	  spouses[spouses[aPeripheral][0]].Delete(idx);
	  spouses[aPeripheral].Delete(0);
          
	  if(isFinal(peelTo.first))
	    {
	      if(peeled.size()!=famSize-1)
		error("Inconsistency observed!\n");
	      break;
	    }
	  
	  if(isLeaf(spouses[aPeripheral][0]))
	    leaf.Push(spouses[aPeripheral][0]);
	  else if(isPeripheral(spouses[aPeripheral][0]))
	    peripheral.Push(spouses[aPeripheral][0]);
	  else if(isRoof(spouses[aPeripheral][0]))
	    UpdateRoof(roof, spouses[aPeripheral][0]);
	  
	  continue;
	}
      
      if(leaf.Length()>0)
	{
	  aLeaf = leaf[0];
	  leaf.Delete(0);
          
	  peeled.push_back(aLeaf);
          
	  peelFrom.first = aLeaf;
	  peelFrom.second = -1;
	  from.push_back(peelFrom);
	  peelingType.push_back(1);
          
	  peelTo.first = parents[aLeaf][0]; 
	  peelTo.second = parents[aLeaf][1];
	  to.push_back(peelTo);
          
	  int idx = offspring[peelTo.first].Find(aLeaf);
	  if(idx<0) 
	    error("Peeling error for person %s in family %s! Check pedigree structure!!\n", ped->persons[family->path[aLeaf]]->pid.c_str(), family->famid.c_str());
	  offspring[peelTo.first].Delete(idx);
          
	  idx = offspring[peelTo.second].Find(aLeaf);
	  if(idx<0) error("Peeling leaf error: %d!\n", aLeaf);
	  offspring[peelTo.second].Delete(idx);
	  parents[aLeaf].Delete(0);
	  parents[aLeaf].Delete(0);
          
	  if(isRoof(peelTo.first))
	    UpdateRoof(roof, peelTo.first);
	  else
	    {
	      if(isPeripheral(peelTo.first)) 
		{ 
		  peripheral.Push(peelTo.first); 
		}
	      if(isPeripheral(peelTo.second))
		{
		  peripheral.Push(peelTo.second);
		}
	    }
	  
	  int pos = find_element(roof, peelTo);
	  if(pos>0)
	    error("Peeling leaf error!\n");
	  //remove_element(roof, pos);
	}
    }        
}

bool ES_Peeling::isFinal(int index)
{
  if(parents[index].Length()==0 && spouses[index].Length()==0 && offspring[index].Length()==0)
    return(true);
  return(false);
}

bool ES_Peeling::isLeaf(int index)
{
  if(offspring[index].Length()==0 && spouses[index].Length()==0)
    return(true);
  return(false);                    
}

bool ES_Peeling::isPeripheral(int index)
{
  if(offspring[index].Length()==0 && parents[index].Length()==0 && spouses[index].Length()==1)
    return(true);
  return(false);
}

bool ES_Peeling::isRoof(int i)
{
  if(spouses[i].Length()==1 && spouses[spouses[i][0]].Length()==1 && 
     parents[i].Length()==0 && parents[spouses[i][0]].Length()==0 && 
     offspring[i].Length()==1 && offspring[spouses[i][0]].Length()==1)
    return(true);
  return(false);   
}

bool ES_Peeling::UpdateRoof(std::vector<pair<int, int> > &roof, int index)
{
  std::pair<int, int> couple;
  couple.first = index;
  couple.second = spouses[index][0];
  if(find_element(roof, couple)>=0) return(false);
  roof.push_back(couple);
  return(true);
}

int ES_Peeling::find_element(std::vector<pair<int, int> >&v, std::pair<int, int> &p)
{
  int position = 0;
  std::vector<pair<int, int> >::iterator it;
  for(it=v.begin(); it!=v.end(); it++)
    {
      if(it->first==p.first && it->second==p.second || it->first==p.second && it->second==p.first)
	return(position);
      position++;
    }
  return(-1);
}

bool ES_Peeling::remove_element(std::vector<pair<int, int> > &v, int index)
{
  if(index<0 || index>=v.size()) return(false);
  v.erase(v.begin()+index);
  return(true);
}

void ES_Peeling::PrintPeelingOrder()
{
  for(int i=0; i<from.size(); i++)
    {
      printf("%s", (const char *) ped->persons[family->path[from[i].first]]->pid);
      if(from[i].second>=0)
	printf(",%s", (const char *) ped->persons[family->path[from[i].second]]->pid);
      printf(" --> %s", (const char *) ped->persons[family->path[to[i].first]]->pid);
      if(to[i].second>=0)
	printf(",%s", (const char *) ped->persons[family->path[to[i].second]]->pid);
      printf("\n");
    }
}

String ES_Peeling::GetPID(int i)
{
  if(i<0) return(String("-"));
  return(ped->persons[family->path[i]]->pid);
}


FamilyLikelihoodES::FamilyLikelihoodES()
{
  InitValues();
}

void FamilyLikelihoodES::InitValues()
{
  allele1 = allele2 = -1;
  famIdx = -1;
  famSize = 0;
  ped = NULL;
  family = NULL;
  transmission = NULL;
  transmission_denovo = NULL;
  transmission_BA = NULL;
  aM = NULL;
  gM = NULL;
  genoIdx.resize(3);
}

FamilyLikelihoodES::~FamilyLikelihoodES()
{
if(transmission!=NULL) FreeTransmissionMatrix();
if(transmission_denovo!=NULL) FreeTransmissionMatrix_denovo();
if(transmission_BA!=NULL) FreeTransmissionMatrix_BA();
}

void FamilyLikelihoodES::SetPedigree(Pedigree *pedptr)
{
  ped = pedptr;
}

void FamilyLikelihoodES::SetFamily(Family *famptr)
{
  family = famptr;
  famSize = family->count;
  nFounders = family->founders;
  penetrances.Dimension(famSize, 10);
  penetrances.Zero();
}

void FamilyLikelihoodES::PreparePeeling()
{
  es.SetFamily(family);
  es.SetupConnections();
  es.BuildInitialPeelable();
  //es.PrintPeelable();
  es.BuildPeelingOrder();
  //es.PrintPeelingOrder();
}

void FamilyLikelihoodES::SetFamilyIndex(int idx)
{
  famIdx=idx;
}

bool FamilyLikelihoodES::isPhenotyped(int who)
{ 
  for(int i=1; i<10; i++)
    if(penetrances[who][i] != penetrances[who][0])
      return(true);
  return(false);
}

void FamilyLikelihoodES::SetAlleles(int a1, int a2)
{ 
  allele1 = a1; allele2=a2;
  genoIdx[0] = glfHandler::GenotypeIndex(allele1, allele1);
  genoIdx[1] = glfHandler::GenotypeIndex(allele1, allele2);
  genoIdx[2] = glfHandler::GenotypeIndex(allele2, allele2);
}

void FamilyLikelihoodES::InitializeStates(double freq){}

void FamilyLikelihoodES::SetFounderPriors(double freq)
{
  if(allele1>4 || allele1<1 || allele2>4 || allele2<1) 
    error("Alleles are not set: %d %d\n", allele1, allele2);
  priors.resize(10);
  for(int i=0; i<priors.size(); i++) priors[i] = 0.0;  
  priors[genoIdx[0]]  = freq*freq;
  priors[genoIdx[1]] = 2*freq*(1-freq);
  priors[genoIdx[2]] = (1-freq)*(1-freq);
}

void FamilyLikelihoodES::SetFounderPriors_BA(double freq)
{
  if(allele1>4 || allele1<1 || allele2>4 || allele2<1) 
    error("Alleles are not set: %d %d\n", allele1, allele2);
  priors.resize(3);
  priors[0]  = freq*freq;
  priors[1] = 2*freq*(1-freq);
  priors[2] = (1-freq)*(1-freq);  
}

void FamilyLikelihoodES::SetGenotypeMutationModel(GenotypeMutationModel* genoMutModel)
{
  gM = genoMutModel;
}

void FamilyLikelihoodES::FreeTransmissionMatrix_BA()
{
  for(int i=0; i<3; i++)
    delete[] transmission_BA[i];
  delete[] transmission_BA;
  transmission_BA = NULL;
}

void FamilyLikelihoodES::FreeTransmissionMatrix()
{
  for(int i=0; i<10; i++)
    delete[] transmission[i];
  delete[] transmission;
  transmission = NULL;
}

void FamilyLikelihoodES::FreeTransmissionMatrix_denovo()
{
  if(transmission!=NULL) FreeTransmissionMatrix();
  for(int i=0; i<10; i++)
    delete[] transmission_denovo[i];
  delete[] transmission_denovo;
  transmission_denovo = NULL;
}

void FamilyLikelihoodES::SetTransmissionMatrix()
{
  if(transmission!=NULL) FreeTransmissionMatrix();
  transmission = new vector<double>*[10]; 
  for(int i=0; i<10; i++)
    transmission[i] = new vector<double>[10];
  
  for(int i=0; i<10; i++)
    for(int j=0; j<10; j++)
      transmission[i][j].resize(10);
  
  int idx1, idx2;
  int geno[4];
  
  for(int i=1; i<=4; i++)
    for(int j=i; j<=4; j++)
      {
	idx1 = glfHandler::GenotypeIndex(i, j);
	
	for(int k=1; k<=4; k++)
	  for(int m=k; m<=4; m++)
	    {
	      idx2 = glfHandler::GenotypeIndex(k,m);
	      
	      geno[0] = glfHandler::GenotypeIndex(i,k);
	      geno[1] = glfHandler::GenotypeIndex(i,m);   
	      geno[2] = glfHandler::GenotypeIndex(j,k);  
	      geno[3] = glfHandler::GenotypeIndex(j,m); 
	      
	      for(int t=0; t<4; t++)
		transmission[idx1][idx2][geno[t]] += 0.25;
	    }
      }  
}

void FamilyLikelihoodES::SetTransmissionMatrix_denovo()
{
 if(transmission_denovo!=NULL) FreeTransmissionMatrix_denovo();
 if(transmission==NULL) SetTransmissionMatrix();
 if(gM==NULL) error("Genotype mutation model has not been set yet!\n");

 transmission_denovo = new vector<double>*[10]; 
 for(int i=0; i<10; i++)
  transmission_denovo[i] = new vector<double>[10];

 for(int i=0; i<10; i++)
  for(int j=0; j<10; j++)
   transmission_denovo[i][j].resize(10);

 for(int i=0; i<10; i++)
  for(int j=0; j<10; j++)
   for(int k=0; k<10; k++)
   {
    double sum = .0;
    for(int m=0; m<10; m++)
     sum += transmission[i][j][m] * gM->genoMutMatrix[m][k];
    transmission_denovo[i][j][k] = sum;
   }
}

void FamilyLikelihoodES::SetTransmissionMatrix_BA()
{
  if(transmission_BA!=NULL) FreeTransmissionMatrix_BA();
  transmission_BA = new vector<double>*[3]; 
  for(int i=0; i<3; i++)
    transmission_BA[i] = new vector<double>[3];
  
  
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      transmission_BA[i][j].resize(3);
  
  transmission_BA[0][0][0] = 1.0;    transmission_BA[0][0][1] = 0.0;    transmission_BA[0][0][2] = 0.0;
  transmission_BA[0][1][0] = 0.5;    transmission_BA[0][1][1] = 0.5;    transmission_BA[0][1][2] = 0.0;
  transmission_BA[0][2][0] = 0.0;    transmission_BA[0][2][1] = 1.0;    transmission_BA[0][2][2] = 0.0;
  transmission_BA[1][0][0] = 0.5;    transmission_BA[1][0][1] = 0.5;    transmission_BA[1][0][2] = 0.0;   
  transmission_BA[1][1][0] = .25;    transmission_BA[1][1][1] = 0.5;    transmission_BA[1][1][2] = .25;   
  transmission_BA[1][2][0] = 0.0;    transmission_BA[1][2][1] = 0.5;    transmission_BA[1][2][2] = 0.5; 
  transmission_BA[2][0][0] = 0.0;    transmission_BA[2][0][1] = 1.0;    transmission_BA[2][0][2] = 0.0;   
  transmission_BA[2][1][0] = 0.0;    transmission_BA[2][1][1] = 0.5;    transmission_BA[2][1][2] = 0.5;
  transmission_BA[2][2][0] = 0.0;    transmission_BA[2][2][1] = 0.0;    transmission_BA[2][2][2] = 1.0;
}

void FamilyLikelihoodES::PrintTransmissionMatrix()
{
  for(int i=0; i<10; i++)
    {
      for(int j=0; j<10; j++)
	{
	  for(int k=0; k<10; k++)
	    printf(" %f", transmission[i][j][k]);
	  printf("\n");
	}
      printf("\n");
    } 
}

void FamilyLikelihoodES::PrintTransmissionMatrix_BA()
{
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<3; j++)
	{
	  for(int k=0; k<3; k++)
	    printf(" %f", transmission_BA[i][j][k]);
	  printf("\n");
	}
      printf("\n");
    } 
}

void FamilyLikelihoodES::FillPenetrance(Family *){}

void FamilyLikelihoodES::FillPenetrance(){}

void FamilyLikelihoodES::PrintPenetrance()
{
  for(int i=0; i<famSize; i++)
    {
      for(int j=0; j<10; j++)
	printf(" %f", penetrances[i][j]);
      printf("\n");
    }
}


void FamilyLikelihoodES::PrintPenetrance_BA()
{
  for(int i=0; i<famSize; i++)
    {
      for(int j=0; j<3; j++)
        printf(" %f", penetrances[i][genoIdx[j]]);
      printf("\n");
    }
}

void FamilyLikelihoodES::FillPenetrance(int){}

double FamilyLikelihoodES::CalculateLikelihood()
{
  marriage_partials.clear();


  for(int i=0; i<es.from.size(); i++)
    switch(es.peelingType[i])
      { 
      case 1:  { peelOffspring2Parents(i); break; }
      case 2:  { peelSpouse2Spouse(i); break; }
      case 3:  { peelParents2Offspring(i); break; }
      default: { error("Peeling type error!\n"); }
      }
  
  int final = es.to[es.from.size()-1].first;
  double lk = 0.0;
  
  for(int i=0; i<10; i++)
    lk+=partials[final][i];
  
  return(lk);
}


double FamilyLikelihoodES::CalculateLikelihood_BA()
{
  marriage_partials.clear();
  
  for(int i=0; i<es.from.size(); i++)
    switch(es.peelingType[i])
      { 
      case 1:  { peelOffspring2Parents_BA(i); break; }
      case 2:  { peelSpouse2Spouse_BA(i); break; }
      case 3:  { peelParents2Offspring_BA(i); break; }
      default: { error("Peeling type error!\n"); }
      }
  
  int final = es.to[es.from.size()-1].first;
  double lk = 0.0;
  
  for(int i=0; i<3; i++)
    lk+=partials[final][i];
  
  return(lk);
}

// family likelihood allowing for de novo mutations
double FamilyLikelihoodES::CalculateLikelihood_denovo()
{
  if(gM==NULL) error("Genotype mutation model is not set!\n");

  marriage_partials.clear();
  
  for(int i=0; i<es.from.size(); i++)
    switch(es.peelingType[i])
      { 
      case 1:  { peelOffspring2Parents_denovo(i); break; }
      case 2:  { peelSpouse2Spouse_denovo(i); break; }
      case 3:  { peelParents2Offspring_denovo(i); break; }
      default: { error("Peeling type error!\n"); }
      }
  
  int final = es.to[es.from.size()-1].first;
  double lk = 0.0;
  
  for(int i=0; i<10; i++)
    lk+=partials[final][i];
  
  return(lk);
}

void FamilyLikelihoodES::peelOffspring2Parents(int idx)
{
  int offspring = es.from[idx].first;
  int fa = es.to[idx].first;
  int mo = es.to[idx].second;
  
  std::map<pair<int, int>, vector<vector<double> > >::iterator it;
  it = marriage_partials.find(es.to[idx]);
  if(it==marriage_partials.end())
    SetMarriagePartials(es.to[idx]);
  
  double partial_lk_sum = 0;
  
  for(int i=0; i<10; i++)
    for(int j=0; j<10; j++)
      {
	partial_lk_sum = 0;
	for(int k=0; k<10; k++)
	  {
	    partial_lk_sum += transmission[i][j][k]*partials[offspring][k];
	  }
	marriage_partials[es.to[idx]][i][j] *= partial_lk_sum;
      }
}


void FamilyLikelihoodES::peelOffspring2Parents_BA(int idx)
{
  int offspring = es.from[idx].first;
  int fa = es.to[idx].first;
  int mo = es.to[idx].second;
  
  std::map<pair<int, int>, vector<vector<double> > >::iterator it;
  it = marriage_partials.find(es.to[idx]);
  if(it==marriage_partials.end())
    SetMarriagePartials(es.to[idx]);
  
  double partial_lk_sum = 0;
  
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      {
	partial_lk_sum = 0;
	for(int k=0; k<3; k++)
	  {
	    partial_lk_sum += transmission_BA[i][j][k]*partials[offspring][k];
	  }
	marriage_partials[es.to[idx]][i][j] *= partial_lk_sum;
      }
}

void FamilyLikelihoodES::peelSpouse2Spouse(int idx)
{
  int spouse_from = es.from[idx].first;
  int spouse_to   = es.to[idx].first;
  
  pair<int, int> spouse;
  bool fa2mo = true;
  if(ped->persons[family->path[spouse_from]]->sex==2)
    { 
      spouse.first = spouse_to;
      spouse.second = spouse_from;
      fa2mo = false;
    }
  else {
    spouse.first = spouse_from;
    spouse.second = spouse_to;
    fa2mo = true;
  }
  
  double partial_lk_sum = 0.0;
  if(marriage_partials.find(spouse)==marriage_partials.end())
    {
      for(int i=0; i<10; i++)
	{
	  partial_lk_sum = 0.0;
	  for(int j=0; j<10; j++)
	    partial_lk_sum += partials[spouse_from][j] ;
	  partials[spouse_to][i] *= partial_lk_sum;
	}
    }
  else if(fa2mo==true){
    for(int i=0; i<10; i++)
      {
	partial_lk_sum = 0.0;
	for(int j=0; j<10; j++)
	  partial_lk_sum += partials[spouse_from][j] * marriage_partials[spouse][j][i] ;
	partials[spouse_to][i] *= partial_lk_sum;
      }
  }
  else {
    for(int i=0; i<10; i++)
      {
	partial_lk_sum = 0.0;
	for(int j=0; j<10; j++)
	  partial_lk_sum += partials[spouse_from][j] * marriage_partials[spouse][i][j] ;
	partials[spouse_to][i] *= partial_lk_sum;
      }
  }
}

void FamilyLikelihoodES::peelSpouse2Spouse_BA(int idx)
{
  int spouse_from = es.from[idx].first;
  int spouse_to   = es.to[idx].first;
  
  pair<int, int> spouse;
  bool fa2mo = true;
  if(ped->persons[family->path[spouse_from]]->sex==2)
    { 
      spouse.first = spouse_to;
      spouse.second = spouse_from;
      fa2mo = false;
    }
  else {
    spouse.first = spouse_from;
    spouse.second = spouse_to;
    fa2mo = true;
  }
  
  double partial_lk_sum = 0.0;
  if(marriage_partials.find(spouse)==marriage_partials.end())
    {
      for(int i=0; i<3; i++)
	{
	  partial_lk_sum = 0.0;
	  for(int j=0; j<3; j++)
	    partial_lk_sum += partials[spouse_from][j] ;
	  partials[spouse_to][i] *= partial_lk_sum;
	}
    }
  else if(fa2mo==true){
    for(int i=0; i<3; i++)
      {
	partial_lk_sum = 0.0;
	for(int j=0; j<3; j++)
	  partial_lk_sum += partials[spouse_from][j] * marriage_partials[spouse][j][i] ;
	partials[spouse_to][i] *= partial_lk_sum;
      }
  }
  else {
    for(int i=0; i<3; i++)
      {
	partial_lk_sum = 0.0;
	for(int j=0; j<3; j++)
	  partial_lk_sum += partials[spouse_from][j] * marriage_partials[spouse][i][j] ;
	partials[spouse_to][i] *= partial_lk_sum;
      }
  }
}

void FamilyLikelihoodES::peelParents2Offspring(int idx)
{
  int fa = es.from[idx].first;
  int mo = es.from[idx].second;
  int offspring = es.to[idx].first;
  
  double partial_lk_sum = 0.0;
  map<pair<int, int>, vector<vector<double> > >::iterator it;
  it = marriage_partials.find(es.from[idx]);
  
  for(int k=0; k<10; k++)
    {
      partial_lk_sum = 0.0;
      
      if(it==marriage_partials.end())
	for(int i=0; i<10; i++)
	  for(int j=0; j<10; j++)
	    partial_lk_sum += partials[fa][i] * partials[mo][j] * transmission[i][j][k];
      
      else
	for(int i=0; i<10; i++)
	  for(int j=0; j<10; j++)
	    partial_lk_sum += partials[fa][i] * marriage_partials[es.from[idx]][i][j] * partials[mo][j] * transmission[i][j][k];  
      
      partials[offspring][k] *= partial_lk_sum;
      
    } 
}

void FamilyLikelihoodES::peelParents2Offspring_BA(int idx)
{
  int fa = es.from[idx].first;
  int mo = es.from[idx].second;
  int offspring = es.to[idx].first;
  
  double partial_lk_sum = 0.0;
  map<pair<int, int>, vector<vector<double> > >::iterator it;
  it = marriage_partials.find(es.from[idx]);
  
  for(int k=0; k<3; k++)
    {
      partial_lk_sum = 0.0;
      
      if(it==marriage_partials.end())
	for(int i=0; i<3; i++)
	  for(int j=0; j<3; j++)
	    partial_lk_sum += partials[fa][i] * partials[mo][j] * transmission_BA[i][j][k];
      
      else
	for(int i=0; i<3; i++)
	  for(int j=0; j<3; j++)
	    partial_lk_sum += partials[fa][i] * marriage_partials[es.from[idx]][i][j] * partials[mo][j] * transmission_BA[i][j][k];  
      
      partials[offspring][k] *= partial_lk_sum;
      
    } 
}

// allowing for de novo mutations while peeling
void FamilyLikelihoodES::peelOffspring2Parents_denovo(int idx)
{
  int offspring = es.from[idx].first;
  int fa = es.to[idx].first;
  int mo = es.to[idx].second;
  
  std::map<pair<int, int>, vector<vector<double> > >::iterator it;
  it = marriage_partials.find(es.to[idx]);
  if(it==marriage_partials.end())
    SetMarriagePartials(es.to[idx]);
  
  double partial_lk_sum = 0;
  
  for(int i=0; i<10; i++)
    for(int j=0; j<10; j++)
      {
	partial_lk_sum = 0;
	for(int k=0; k<10; k++)
	  //for(int m=0; m<10; m++)
	    //partial_lk_sum += transmission[i][j][m] * gM->genoMutMatrix[m][k] * partials[offspring][k];
	  partial_lk_sum += transmission_denovo[i][j][k] * partials[offspring][k]; //pre-calculated the de novo transition matrix
	marriage_partials[es.to[idx]][i][j] *= partial_lk_sum;
      }
}

void FamilyLikelihoodES::peelSpouse2Spouse_denovo(int idx)
{
  int spouse_from = es.from[idx].first;
  int spouse_to   = es.to[idx].first;
  
  pair<int, int> spouse;
  bool fa2mo = true;
  if(ped->persons[family->path[spouse_from]]->sex==2)
    { 
      spouse.first = spouse_to;
      spouse.second = spouse_from;
      fa2mo = false;
    }
  else {
    spouse.first = spouse_from;
    spouse.second = spouse_to;
    fa2mo = true;
  }
  
  double partial_lk_sum = 0.0;
  if(marriage_partials.find(spouse)==marriage_partials.end())
    {
      for(int i=0; i<10; i++)
	{
	  partial_lk_sum = 0.0;
	  for(int j=0; j<10; j++)
	    partial_lk_sum += partials[spouse_from][j];
	  partials[spouse_to][i] *= partial_lk_sum;
	}
    }
  else if(fa2mo==true){
    for(int i=0; i<10; i++)
      {
	partial_lk_sum = 0.0;
	for(int j=0; j<10; j++)
	  partial_lk_sum += partials[spouse_from][j] * marriage_partials[spouse][j][i] ;
	partials[spouse_to][i] *= partial_lk_sum;
      }
  }
  else {
    for(int i=0; i<10; i++)
      {
	partial_lk_sum = 0.0;
	for(int j=0; j<10; j++)
	  partial_lk_sum += partials[spouse_from][j] * marriage_partials[spouse][i][j] ;
	partials[spouse_to][i] *= partial_lk_sum;
      }
  }
}

void FamilyLikelihoodES::peelParents2Offspring_denovo(int idx)
{
  int fa = es.from[idx].first;
  int mo = es.from[idx].second;
  int offspring = es.to[idx].first;
  
  double partial_lk_sum = 0.0;
  map<pair<int, int>, vector<vector<double> > >::iterator it;
  it = marriage_partials.find(es.from[idx]);
  
  for(int k=0; k<10; k++)
    {
      partial_lk_sum = 0.0;
      if(it==marriage_partials.end())
	{
	  for(int i=0; i<10; i++)
	    for(int j=0; j<10; j++)
	      //for(int m=0; m<10; m++)
		//partial_lk_sum += partials[fa][i] * partials[mo][j] * transmission[i][j][m] * gM->genoMutMatrix[m][k];
		partial_lk_sum += partials[fa][i] * partials[mo][j] * transmission_denovo[i][j][k];
	}
      else
	{
	  for(int i=0; i<10; i++)
	    for(int j=0; j<10; j++)
	      //for(int m=0; m<10; m++)
		//partial_lk_sum += partials[fa][i] * marriage_partials[es.from[idx]][i][j] * partials[mo][j] * transmission[i][j][m] * gM->genoMutMatrix[m][k];  
		partial_lk_sum += partials[fa][i] * marriage_partials[es.from[idx]][i][j] * partials[mo][j] * transmission[i][j][k];  
	}
      partials[offspring][k] *= partial_lk_sum;
    } 
}

void FamilyLikelihoodES::CalculateLikelihood(Family *){}
void FamilyLikelihoodES::CalculateLikelihood(int){}

void FamilyLikelihoodES::SetMarriagePartials(std::pair<int, int>& parents)
{
  int fa = parents.first;
  int mo = parents.second;
  
  vector<vector<double> > partial;
  partial.resize(10);
  for(int i=0; i<10; i++)
    partial[i].resize(10);
  
  for(int i=0; i<10; i++)
    for(int j=0; j<10; j++)
      partial[i][j] = 1.0;
  
  marriage_partials[parents] = partial;    
}

void FamilyLikelihoodES::SetMarriagePartials_BA(std::pair<int, int>& parents)
{
  int fa = parents.first;
  int mo = parents.second;
  
  vector<vector<double> > partial;
  partial.resize(3);
  for(int i=0; i<3; i++)
    partial[i].resize(3);
  
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      partial[i][j] = 1.0;
  
  marriage_partials[parents] = partial;    
}

void FamilyLikelihoodES::InitializePartials()
{
  partials.Dimension(famSize, 10);  
  partials.Zero();  
  for(int i=0; i<famSize; i++)
    {
     if(ped->persons[family->path[i]]->isFounder())
      for(int j=0; j<10; j++)
        partials[i][j] = priors[j]*penetrances[i][j];
     else
      for(int j=0; j<10; j++)
	  partials[i][j] = penetrances[i][j];
    }   
}


void FamilyLikelihoodES::InitializePartials_BA()
{
  partials.Dimension(famSize, 3);  
  partials.Zero();
  for(int i=0; i<famSize; i++)
    {
      if(ped->persons[family->path[i]]->isFounder())
        for(int j=0; j<3; j++)
         partials[i][j] = priors[j]*penetrances[i][genoIdx[j]];
      else
        for(int j=0; j<3; j++)
	  partials[i][j] = penetrances[i][genoIdx[j]];
    }
}

String FamilyLikelihoodES::GetPID(int i)
{
  if(i<0) return(String("-"));
  return(ped->persons[family->path[i]]->pid);
}
