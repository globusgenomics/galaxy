#include "NucFamGenotypeLikelihood.h"
#include "StringMap.h"
#include <stdlib.h>
#include <iostream>
#include "CmdLinePar.h"
#include "MutationModel.h"
#include <map>
#include "FamilyLikelihoodES.h"
#include "FamilyLikelihoodSeq.h"

StringArray tokens;
StringArray lists;

int readGLFannoFile(String GLFannoFile, StringMap *glfMap)
{
  IFILE fh = ifopen(GLFannoFile.c_str(), "r");
  if(fh==NULL)
    error("%s open failed\n", GLFannoFile.c_str());
  
  String line;
  StringArray tokens;
  
  int count = 0;
  while(!ifeof(fh))
    {
      line.ReadLine(fh);
      tokens.Clear();
      tokens.ReplaceTokens(line);
      if(tokens.Length()<2) { continue; }
      lists.Add(tokens[1]);
      glfMap->Add(tokens[0],&lists[count]);
      count++;
    }
  ifclose(fh);
  return(count);
}

int LoadPositionFile(String &file, std::map<String, int>& positionMap)
{
  FILE *fh = fopen(file.c_str(), "r");
  if(fh==NULL) error("Open position file %s failed!\n", file.c_str());
  String buffer;
  StringArray tokens;
  while(!feof(fh))
    {
      buffer.ReadLine(fh);
      tokens.ReplaceTokens(buffer);
      if(tokens.Length()==0) continue;
      String chrPos = tokens[0]+":"+tokens[1];
      positionMap[chrPos]++;     
    }
  fclose(fh);
  return(0);
}

int main(int argc, char * argv[])
{
  double posterior = 0.5;
  int minTotalDepth = 0;
  int maxTotalDepth = 0;
  double minPS = 0;
  int minMapQuality = 0;
  String pedFile, datFile, glfListFile;
  String vcfFile = "variantCalls.vcf";
  String positionfile;
  String chrs2process;
  double theta = 0.001;
  double tstv_ratio = 2.0;
  double precision = 0.0001;
  int num_threads = 1;
  bool denovo = false;
  double denovo_mut_rate = 1.5e-08;
  double denovo_tstv_ratio = 2.0;
  double denovoLR = 1;
  bool gl_off = false;

  ParameterList pl;

  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Map Quality Filter")
    LONG_INTPARAMETER("minMapQuality", &minMapQuality)
    LONG_PARAMETER_GROUP("Depth Filter")
    LONG_INTPARAMETER("minDepth", &minTotalDepth)
    LONG_INTPARAMETER("maxDepth", &maxTotalDepth)
    LONG_DOUBLEPARAMETER("minPercSampleWithData", &minPS)
    LONG_PARAMETER_GROUP("Scaled mutation rate")
    LONG_DOUBLEPARAMETER("theta", &theta)
    LONG_PARAMETER_GROUP("Prior of ts/tv ratio")
    LONG_DOUBLEPARAMETER("poly_tstv", &tstv_ratio)
  LONG_PARAMETER_GROUP("de novo mutation")
      LONG_PARAMETER("denovo", &denovo)
      LONG_DOUBLEPARAMETER("rate_denovo", &denovo_mut_rate)
      LONG_DOUBLEPARAMETER("tstv_denovo", &denovo_tstv_ratio)
      LONG_DOUBLEPARAMETER("minLLR_denovo", &denovoLR)
    LONG_PARAMETER_GROUP("Optimization precision")
    LONG_DOUBLEPARAMETER("prec", &precision)
  LONG_PARAMETER_GROUP("Multiple threading")
    LONG_INTPARAMETER("nthreads", &num_threads)
  LONG_PARAMETER_GROUP("Chromosomes to process")
    LONG_STRINGPARAMETER("chr2process", &chrs2process)
  LONG_PARAMETER_GROUP("Output")
    LONG_STRINGPARAMETER("vcf", &vcfFile)
    LONG_PARAMETER("gl_off", &gl_off)
    END_LONG_PARAMETERS();
  
  pl.Add(new StringParameter('p', "pedfile", pedFile));
  pl.Add(new StringParameter('d', "datfile", datFile));
  pl.Add(new StringParameter('g', "glfIndexFile", glfListFile));
  pl.Add(new DoubleParameter('c', "posterior cutoff", posterior));
  pl.Add(new LongParameters("Additional Options", longParameters));
  pl.Read(argc, argv);

  pl.Status();


  if(pedFile.Length()==0)
    error("pedFile not provided for input!\n");
  if(datFile.Length()==0)
    error("datFile not provided for input!\n");
  if(glfListFile.Length()==0)
    error("glfListFile not provided for input!\n");
  if(vcfFile.Length()==0)
    error("vcfFile not provided for output!\n");

  std::map<String, int> positionMap;
  if(positionfile.Length()>0) LoadPositionFile(positionfile, positionMap);

  #ifdef _OPENMP
  if(num_threads>0) omp_set_num_threads(num_threads);
  #endif

  CmdLinePar par;
  par.minTotalDepth = minTotalDepth;
  par.maxTotalDepth = maxTotalDepth;
  par.minMapQuality = minMapQuality;
  par.minPS = minPS;
  par.posterior     = posterior;
  par.precision     = precision;
  par.denovo_mut_rate = denovo_mut_rate;
  par.denovo_tstv_ratio = denovo_tstv_ratio;
  par.denovo = denovo;
  par.denovoLR = denovoLR;
  par.gl_off = gl_off;
  par.chrs2process = chrs2process;

  if(denovo && denovoLR<0) error("denovo_min_LLR can only be greater than 0 !\n");

  // Prior of ts and tv
  double prior_ts = tstv_ratio/(tstv_ratio + 1);
  double prior_tv = (1-prior_ts)/2;

  double polyPrior;
  double refTransition = .0;
  double refTransvers1 = .0;
  double refTransvers2 = .0;
  double tstvs1=0; double tstvs2=.0; double tvs1tvs2=.0;
  double refTransitionFreq, refTransvers1Freq, refTransvers2Freq, tstvs1Freq, tstvs2Freq, tvs1tvs2Freq;

  StringMap glfMap;
  String glfFileKey;
  readGLFannoFile(glfListFile, &glfMap);

  Pedigree ped;
  PedigreeGLF pedGLF;

  IFILE datFH = ifopen(datFile, "r");
  IFILE pedFH = ifopen(pedFile, "r");
  FILE *vcfFH = fopen(vcfFile, "w");
  if(datFH==NULL)
    error("datFile open for input failed!\n");
  if(pedFH==NULL)
    error("pedFile open for input failed!\n");
  if(vcfFH==NULL)
    error("vcfFile can not be opened for output!\n");

  ped.Prepare(datFH);
  ped.Load(pedFH);

  if(datFH != NULL) ifclose(datFH);
  if(pedFH != NULL) ifclose(pedFH);

//  This is for debug purpose which is to compare results of nuclear families
//      using two algorithms. They should give IDENTICAL results
//  for(int f=0; f<ped.familyCount; f++)
//   ped.families[f]->generations=3;


  pedGLF.SetGLFMap(&glfMap);

  pedGLF.SetPedGLF(&ped);

  FamilyLikelihoodSeq famlk[7];
  for(int i=0; i<7; i++)
    {
      famlk[i].SetCmdLinePar(&par);
      if(denovo) famlk[i].SetDenovoMutationModel();
      famlk[i].SetTheta(theta);
      famlk[i].SetGLF(&pedGLF);
      famlk[i].InitFamilyLikelihoodES();
    }
  
  int cnt=0;
  int cntSec=0;
  int totalEntryCnt = 0;

  //sumarry statistics
  uint minTotalDepthFilter = 0;
  uint maxTotalDepthFilter = 0;
  uint maxAvgDepthFilter = 0;
  uint minAvgDepthFilter = 0;
  uint minMapQualFilter = 0;
  uint minAvgMapQualFilter = 0;
  uint minPSFilter = 0;
  
  int refBaseCounts[5] = {0,0,0,0,0};
  int homoRef = 0;
  int transitions = 0;
  int transversions = 0;
  int otherPolymorphism = 0;
  int tstvs1Cnt = 0;
  int tstvs2Cnt = 0;
  int tvs1tvs2Cnt = 0;
  int nocall = 0;
  int actuaBases = 0;

  time_t t; time(&t);

  std::map<String, int> chrs2process_map;
  StringArray chrs;
  chrs.AddTokens(chrs2process, ',');
  for(int cidx=0; cidx<chrs.Length(); cidx++)
   chrs2process_map[chrs[cidx]]++;
  
  printf("Analysis started on %s\n", ctime(&t));
  Matrix pen; 
  int maxidx = 0;
  while(pedGLF.Move2NextSection())
    {
     if(chrs2process_map.size()>0 && chrs2process_map[pedGLF.GetNonNULLglf()->label]<1) continue;

      homoRef = transitions = transversions = otherPolymorphism = tstvs1Cnt = tstvs2Cnt = tvs1tvs2Cnt = nocall = actuaBases = 0; for(int k=0; k<5; k++) refBaseCounts[k]=0; 
      cnt=cntSec=totalEntryCnt = 0; minTotalDepthFilter = maxTotalDepthFilter = maxAvgDepthFilter = minAvgDepthFilter = minMapQualFilter = minAvgMapQualFilter = 0;
    
      while(pedGLF.Move2NextBaseEntry())
	{
	  for(int r=0; r<7; r++)
	    famlk[r].FillPenetrance();

	  if(totalEntryCnt==0) totalEntryCnt = pedGLF.GetNonNULLglf()->maxPosition;

	  if(positionfile.Length()>0)
	    {
	      int pos = pedGLF.currentPos+1;
	      String chrPos = pedGLF.GetNonNULLglf()->label+":"+pos;
	      if(positionMap[chrPos]==0) continue;
	    }
	  int refBase  = pedGLF.GetRefBase();
	  if(refBase!=1 && refBase!=2 && refBase!=3 && refBase!=4) continue;
	  refBaseCounts[refBase]++;

	  famlk[0].CalcReadStats();

	  if(famlk[0].totalDepth<minTotalDepth) { minTotalDepthFilter++; continue; }
	  if(maxTotalDepth>0 && famlk[0].totalDepth>maxTotalDepth) { maxTotalDepthFilter++; continue; }
	  if(famlk[0].percSampWithData*100<minPS){ minPSFilter++; continue; }
	  if(famlk[0].avgMapQual<minMapQuality) { minMapQualFilter++; continue; }

	  int ts       = (Poly::ts(refBase));   //transition
	  int tvs1     = (Poly::tvs1(refBase)); //transvertions1
	  int tvs2     = (Poly::tvs2(refBase)); //transvertion2

	  polyPrior = famlk[0].GetPolyPrior();

# ifdef _OPENMP
# pragma omp parallel sections
# endif
	  {
# ifdef _OPENMP
#pragma omp section
# endif
	    {
	      if(par.denovo==false)
	      {
		//Calculate likelihood assuming everyone is homozygous for the reference
		double lRef = log10(1-polyPrior) + famlk[0].MonomorphismLogLikelihood(refBase);
		famlk[0].varllk[0]  = lRef;
		famlk[0].varllk_noprior[0] = lRef - log10(1-polyPrior);
		famlk[0].varfreq[0] = 1.0;
	      }
	      else 
	      {
		//Calculate likelihood assuming everyone is homozygous for the reference allowing for denovo mutations
		double lRef_denovo = log10(1-polyPrior) + famlk[0].MonomorphismLogLikelihood_denovo(refBase, refBase==4 ? refBase-1 : refBase+1);
		famlk[0].varllk[0] = lRef_denovo;
		famlk[0].varllk_noprior[0] = lRef_denovo - log10(1-polyPrior);
		famlk[0].varfreq[0] = 1.0;
	      }
	    }
# ifdef _OPENMP
# pragma omp section
# endif
	    {
	      // Calculate likelihoods for the most likelily SNP configurations
	      refTransition = log10(polyPrior * prior_ts) + famlk[1].PolymorphismLogLikelihood(refBase, ts);
	      refTransitionFreq = famlk[1].GetMinimizer();
	      famlk[0].varllk[1]  = refTransition;
		   famlk[0].varllk_noprior[1] = refTransition -  log10(polyPrior * 2./3. );
	      famlk[0].varfreq[1] = refTransitionFreq;
	    }
# ifdef _OPENMP
# pragma omp section
# endif
	    {
	      refTransvers1 = log10(polyPrior * prior_tv) + famlk[2].PolymorphismLogLikelihood(refBase, tvs1);
	      refTransvers1Freq = famlk[2].GetMinimizer();
	      famlk[0].varllk[2]  = refTransvers1;
	    	famlk[0].varllk_noprior[2] = refTransvers1 - log10(polyPrior * 1./6.);
	      famlk[0].varfreq[2] = refTransvers1Freq;
	    }
# ifdef _OPENMP
# pragma omp section
# endif
	    {
	      refTransvers2 = log10(polyPrior * prior_tv) + famlk[3].PolymorphismLogLikelihood(refBase, tvs2);
	      refTransvers2Freq = famlk[3].GetMinimizer();
	      famlk[0].varllk[3]  = refTransvers2;
	      famlk[0].varllk_noprior[3] = refTransvers2 - log10(polyPrior * 1./6.);
	      famlk[0].varfreq[3] = refTransvers2Freq;
	    }
	  }
	  maxidx = famlk[0].CalcVarPosterior(4);

	  // Calculate likelihoods for less likely SNP configurations
	  if(famlk[0].varPostProb<0.99)
	    {
# ifdef _OPENMP
#pragma omp parallel sections
# endif
	      {
# ifdef _OPENMP
#pragma omp section
# endif
		{
		  tstvs1   = log10(polyPrior * 0.001) + famlk[4].PolymorphismLogLikelihood(ts, tvs1);
		  tstvs1Freq = famlk[4].GetMinimizer();
		  famlk[0].varllk[4]  = tstvs1;
		  famlk[0].varllk_noprior[4] = tstvs1 -  log10(polyPrior * 0.001);
		  famlk[0].varfreq[4] = tstvs1Freq;
		}
# ifdef _OPENMP
#pragma omp section
# endif
		{
		  tstvs2   = log10(polyPrior * 0.001) + famlk[5].PolymorphismLogLikelihood(ts, tvs2);
		  tstvs2Freq = famlk[5].GetMinimizer();
		  famlk[0].varllk[5]  = tstvs2;
		  famlk[0].varllk_noprior[5] = tstvs2 - log10(polyPrior * 0.001) ;
		  famlk[0].varfreq[5] = tstvs2Freq;
		}
# ifdef _OPENMP
#pragma omp section
# endif
		{
		  tvs1tvs2 = log10(polyPrior * 0.001) + famlk[6].PolymorphismLogLikelihood(tvs1, tvs2);
		  tvs1tvs2Freq = famlk[6].GetMinimizer();
		  famlk[0].varllk[6]  = tvs1tvs2;
		  famlk[0].varllk_noprior[6] = tvs1tvs2 - log10(polyPrior * 0.001);
		  famlk[0].varfreq[6] = tvs1tvs2Freq;
		}
	      }
	      maxidx = famlk[0].CalcVarPosterior(7);
	    }

	  if(famlk[0].varPostProb<posterior) { nocall++; continue; }

	  switch(maxidx)
	    {
	    case 0: { homoRef++; famlk[0].SetAlleles(refBase, refBase==4?refBase-1:refBase+1); famlk[0].min=1.0; break;}
	    case 1: { transitions++; famlk[0].SetAlleles(refBase, ts); famlk[0].min=famlk[1].min; break;}
	    case 2: { transversions++; famlk[0].SetAlleles(refBase, tvs1);  famlk[0].min=famlk[2].min; break;}
	    case 3: { transversions++; famlk[0].SetAlleles(refBase, tvs2); famlk[0].min=famlk[3].min; break;}
	    case 4: { tstvs1Cnt++; famlk[0].SetAlleles(ts, tvs1); famlk[0].min=famlk[4].min; break; }
	    case 5: { tstvs2Cnt++; famlk[0].SetAlleles(ts, tvs2); famlk[0].min=famlk[5].min; break; }
	    case 6: { tvs1tvs2Cnt++; famlk[0].SetAlleles(tvs1, tvs2); famlk[0].min=famlk[6].min; break;}
	    case -1: { nocall++; break;}
	    default: error("Invalid maxidx!\n");
	    }

	  if(maxidx==-1 || maxidx==0 && par.denovo==false) continue;

	  if(maxidx==0)
	  {
	    double lk_mono = famlk[0].MonomorphismLogLikelihood(refBase);
	    famlk[0].min = 1.0;
	    famlk[0].denovoLR = famlk[0].varllk_noprior[0] - lk_mono;
	    if(famlk[0].denovoLR <= log10(par.denovoLR)) continue;
	  }
	  else {
	    if(par.denovo)
	    {
	      par.denovo = false;
	      double lk_poly = famlk[0].PolymorphismLogLikelihood(famlk[0].allele1, famlk[0].allele2);
	      famlk[0].denovoLR = famlk[0].varllk_noprior[maxidx] - lk_poly;
	      par.denovo = true;
	    }
	  }

    if(maxidx==0)
    {
      famlk[0].denovo_mono = true;
     	famlk[0].CalcPostProb(1.0);
    }
    else
	   famlk[0].CalcPostProb(famlk[0].min); //mlefreq is the freq of the reference allele

	  denovo ? famlk[0].OutputVCF_denovo(vcfFH) : famlk[0].OutputVCF(vcfFH);
	  famlk[0].denovo_mono = false;
	}
      //output summary statisitcs after the run
      int totalBases = 0;
      for(int i=0; i<5; i++){
	totalBases += refBaseCounts[i];
      }
      printf("Summary of reference -- %s\n", pedGLF.GetNonNULLglf()->label.c_str()); 
      printf("Total Entry Count: %9d\n", totalEntryCnt);
      printf("Total Base Cout: %9d\n", totalBases);
      printf("Total '0' Base Count: %9d\n", refBaseCounts[0]);
      printf("Non-Polymorphic Count: %9d\n", homoRef);
      printf("Transition Count: %9d\n", transitions);
      printf("Transversion Count: %9d\n", transversions);
      printf("Other Polymorphism Count: %9d\n", tstvs1Cnt+tstvs2Cnt+tvs1tvs2Cnt);
      printf("Filter counts:\n");
      printf("\tminMapQual %u\n", minMapQualFilter);
      printf("\tminTotalDepth %u\n", minTotalDepthFilter);
      printf("\tmaxTotalDepth %u\n", maxTotalDepthFilter);
      printf("Hard to call: %9d\n", nocall);
      printf("Skipped bases: %u\n", totalEntryCnt-homoRef-transitions-transversions-(tstvs1Cnt+tstvs2Cnt+tvs1tvs2Cnt));
      time_t tfinished = t;
      time(&tfinished);
      time_t duration = tfinished - t;
      printf("Analysis ended on %s\n", ctime(&tfinished));
      printf("Running time is %u seconds\n\n", (unsigned int)(duration));
    }
  fclose(vcfFH);
  return(0);
}
