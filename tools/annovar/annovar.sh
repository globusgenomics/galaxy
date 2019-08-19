#!/bin/bash

test="N"
dofilter="N"

#########################
#	   DEFINE SOME
#	    FUNCTIONS
#########################

function usage(){
	echo "usage: $0 todo"
}

function runfilter(){
	ifile=$1	
	columnname=$2
	threshold=$3

	if [[ $threshold == "-1" ]]
	then
		echo "not filtering"
		return
	fi
	
	echo "filtering: $columnname, $threshold"
	cat $ifile

	#get column number corresponding to column header
	column=`awk 'BEGIN{
					FS="\t";
					col=-1
				}{
					if(FNR==1){
						for(i=1;i<=NF;i++){
							if($i == "'"${columnname}"'") 
								col=i 
						} 
						print col 
					}
				}' $ifile `

	if [ $column == -1 ]
	then
		echo "no such column, exiting"
		return
	fi	

	#perform filtering using the threshold
	awk 'BEGIN{
		FS="\t";
		OFS="\t";
	}{
		if(FNR==1) 
			print $0; 
		if(FNR>1){
			if( $"'"${column}"'" == "" )  # empty column, then print
				print $0
			else if ("'"${threshold}"'" == "text"){}  #if set to text dont check threshold
				
			else if ($"'"${column}"'" < "'"${threshold}"'")  #else do check it
				print $0	
		}
	}' $ifile > tmpfile

	mv tmpfile $ifile	
}

# arguments: originalfile,resultfile,chrcol,startcol,endcol,refcol,obscol,addcols
function joinresults(){
	ofile=$1
	rfile=$2
	colchr=$3
	colstart=$4
	colend=$5
	colref=$6
	colobs=$7
	addcols=$8 #e.g. "B.col1,B.col2"
	
	test="N"
	
	# echo "joining result with original file"
	if [ $test == "Y" ]
	then 	
		echo "ofile: $ofile"
		head $ofile 
		echo "rfile: $rfile"
		head $rfile
	fi
	numlines=`wc $rfile | cut -d" " -f2`
	
	# if empty results file, just add header fields
	if [[ ! -s $rfile ]] 
	then			
		dummycol=${addcols:2}
		outputcol=${dummycol//",B."/"	"}
		numcommas=`echo "$addcols" | grep -o "," | wc -l`		
		
		awk 'BEGIN{FS="\t";OFS="\t"}{
				if(FNR==1)
					print $0,"'"$outputcol"'"; 
				else{
					printf $0
					for(i=0;i<="'"$numcommas"'"+1;i++)
						printf "\t"
					printf "\n"
				}
			}END{}' $ofile > tempofile
			
			mv tempofile $ofile		
		return
	fi
	

	#get input file column names for cgatools join
	col_chr_name=`head -1  $rfile | cut -f${colchr}`
	col_start_name=`head -1  $rfile | cut -f${colstart}`
	col_end_name=`head -1  $rfile | cut -f${colend}`
	col_ref_name=`head -1  $rfile | cut -f${colref}`
	col_obs_name=`head -1  $rfile | cut -f${colobs}`

	#get annotation file column names for cgatools join
	chr_name=`head -1  $ofile | cut -f${chrcol}`
	start_name=`head -1  $ofile | cut -f${startcol}`
	end_name=`head -1  $ofile | cut -f${endcol}`
	ref_name=`head -1  $ofile | cut -f${refcol}`
	obs_name=`head -1  $ofile | cut -f${obscol}`

	if [ $test == "Y" ]	
	then
		echo "input file"
		echo "chr   col: $col_chr_name ($colchr)"	
		echo "start col: $col_start_name ($colstart)"	
		echo "end   col: $col_end_name ($colend)"	
		echo "ref   col: $col_ref_name ($colref)"	
		echo "obs   col: $col_obs_name ($colobs)"	
		echo ""
		echo "annotation file"
		echo "chr   col: $chr_name ($chrcol)"	
		echo "start col: $start_name ($startcol)"	
		echo "end   col: $end_name ($endcol)"	
		echo "ref   col: $ref_name ($refcol)"	
		echo "obs   col: $obs_name ($obscol)"	
	fi

	#perform join
	cgatools join --beta \
		--input $ofile $rfile \
		--output temporiginal \
		--match ${chr_name}:${col_chr_name} \
		--match ${start_name}:${col_start_name} \
		--match ${end_name}:${col_end_name} \
		--match ${ref_name}:${col_ref_name} \
		--match ${obs_name}:${col_obs_name} \
		--select A.*,$addcols \
		--always-dump \
		--output-mode compact 

	#replace originalfile
	sed -i 's/^>//g' temporiginal #join sometimes adds a '>' symbol to header
	mv temporiginal originalfile
		
	if [ $test == "Y" ]
	then
		echo "joining complete"
		head originalfile
		echo ""	
	fi
	
}




#################################
#
#	   PARSE PARAMETERS
#
#################################


set -- `getopt -n$0 -u -a --longoptions="inputfile: buildver: humandb: varfile: VCF: chrcol: startcol: endcol: refcol: obscol: vartypecol: convertcoords: geneanno: hgvs: verdbsnp: tfbs: mce: cytoband: segdup: dgv: gwas: ver1000g: cg46: cg69: impactscores: newimpactscores: otherinfo: esp: gerp: cosmic65: cosmic67: cosmic70: clinvar: nci60: outall: outfilt: outinvalid: scriptsdir: dorunannovar: dofilter: filt_dbsnp: filt1000GALL: filt1000GAFR: filt1000GAMR: filt1000GASN: filt1000GEUR: filtESP6500ALL: filtESP6500EA: filtESP6500AA: filtcg46: filtcg69: dummy:" "h:" "$@"` || usage
[ $# -eq 0 ] && usage



while [ $# -gt 0 ]
do
    case "$1" in
       	--inputfile)      			infile=$2;shift;;  # inputfile
		--buildver)					buildvertmp=$2;shift;; # hg18 or hg19
		--humandb)					humandbtmp=$2;shift;; # location of humandb database
	 	--varfile)      			varfile=$2;shift;; # Y or N  
		--VCF)						vcf=$2;shift;; #Y or N
		--chrcol)      				chrcol=$2;shift;;  # which column has chr 
		--startcol)      			startcol=$2;shift;;  # which column has start
		--endcol)      				endcol=$2;shift;;  # which column has end
		--refcol)      				refcol=$2;shift;;  # which column has ref
		--obscol)      				obscol=$2;shift;;  # which column has alt
		--vartypecol)      			vartypecol=$2;shift;;  # which column has vartype
		--convertcoords)			convertcoords=$2;shift;;  # Y or N convert coordinate from CG to 1-based?		
		--geneanno)      			geneanno=$2;shift;; # comma-separated list of strings refSeq, knowngene, ensgene  
		--hgvs)						hgvs=$2;shift;;
		--verdbsnp)					verdbsnp=$2;shift;; #comma-separated list of dbsnp version to annotate with (e.g. "132,135NonFlagged,137")"
		--tfbs)      				tfbs=$2;shift;; 	# Y or N 
		--mce)      				mce=$2;shift;; 	# Y or N 
		--cytoband)      			cytoband=$2;shift;; # Y or N  
		--segdup)      				segdup=$2;shift;; 	# Y or N 				
        --dgv)      				dgv=$2;shift;; 	# Y or N 
		--gwas)      				gwas=$2;shift;; 	# Y or N 
		--ver1000g) 				ver1000g=$2;shift;; 	# Y or N 
		--cg46)						cg46=$2;shift;;
		--cg69)						cg69=$2;shift;;		
		--impactscores)      		impactscores=$2;shift;; # Y or N 
		--newimpactscores)      	newimpactscores=$2;shift;; # Y or N 
		--otherinfo)				otherinfo=$2;shift;; 
		--scriptsdir)	      		scriptsdirtmp=$2;shift;; # Y or N 
		--esp)      				esp=$2;shift;; 	# Y or N 
		--gerp)      				gerp=$2;shift;; 	# Y or N 
		--cosmic65)					cosmic65=$2;shift;;  # Y or N
		--cosmic67)					cosmic67=$2;shift;;  # Y or N
                --cosmic70)                                     cosmic70=$2;shift;;  # Y or N
		--nci60)					nci60=$2;shift;;  # Y or N
		--clinvar)					clinvar=$2;shift;;  # Y or N
		--filt_dbsnp)				filt_dbsnp=$2;shift;;
		--filt1000GALL)				threshold_1000g_ALL=$2;shift;; #threshold value
		--filt1000GAFR)				threshold_1000g_AFR=$2;shift;; #threshold value
		--filt1000GAMR)				threshold_1000g_AMR=$2;shift;; #threshold value
		--filt1000GASN)				threshold_1000g_ASN=$2;shift;; #threshold value
		--filt1000GEUR)				threshold_1000g_EUR=$2;shift;; #threshold value
		--filtESP6500ALL)			threshold_ESP6500_ALL=$2;shift;; #threshold value
		--filtESP6500EA)			threshold_ESP6500_EA=$2;shift;; #threshold value
		--filtESP6500AA)			threshold_ESP6500_AA=$2;shift;; #threshold value
		--filtcg46)					threshold_cg46=$2;shift;;
		--filtcg69)					threshold_cg69=$2;shift;;
		--outall)      				outfile_all=$2;shift;; # file 
		--outfilt)      			outfile_filt=$2;shift;; # file
		--outinvalid)				outfile_invalid=$2;shift;; #file
		--dorunannovar)				dorunannovar=$2;shift;; 	#Y or N		
       -h)        	shift;;
	   --)        	shift;break;;
       -*)        	usage;;
       *)         	break;;            
    esac
    shift
done

#sometimes galaxy screws up these variables after updates, if comma-separated list, use only what is before first comma
humandb=${humandbtmp%,*}
buildver=${buildvertmp%,*}
scriptsdir=${scriptsdirtmp%,*}


if [ $test == "Y" ]
then
	echo "dorunannovar: $dorunannovar"
	echo "infile: $infile"
	echo "buildver: $buildver"
	echo "annovardb: $humandb"
	echo "verdbnsp: $verdbsnp"
	echo "geneanno: $geneanno"
	echo "tfbs: $tfbs"
	echo "mce: $mce"
	echo "cytoband: $cytoband"
	echo "segdup: $segdup"
	echo "dgv: $dgv"
	echo "gwas: $gwas"	
	echo "g1000: ${g1000}"
	echo "cg46: ${cg46}"	
	echo "cg69: ${cg69}"
	echo "impactscores: $impactscores"
	echo "impactscores: $newimpactscores"	
	echo "esp: $esp"
	echo "gerp: $gerp"
	echo "cosmic: $cosmic"
	echo "outfile: $outfile_all"
	echo "outinvalid: $outfile_invalid"
	echo "outfiltered: $outfile_filt"
	echo "varfile: $varfile"
	echo "vcf" $vcf
	echo "chrcol: $chrcol"
	echo "startcol: $startcol"
	echo "endcol: $endcol"
	echo "refcol: $refcol"
	echo "obscol: $obscol"
	echo "convertcoords: $convertcoords"
	echo "vartypecol: $vartypecol"
	echo "dofilter: $dofilter"
	echo "threshold_1000g_ALL  : $threshold_1000g_ALL"
	echo "threshold_1000g_AFR  : $threshold_1000g_AFR"
	echo "threshold_1000g_AMR  : $threshold_1000g_AMR"
	echo "threshold_1000g_ASN  : $threshold_1000g_ASN"
	echo "threshold_1000g_EUR  : $threshold_1000g_EUR"
	echo "threshold_ESP6500_ALL: $threshold_ESP6500_ALL"
	echo "threshold_ESP6500_EA : $threshold_ESP6500_EA"
	echo "threshold_ESP6500_AA : $threshold_ESP6500_AA"

fi



############################################
#
#       Annotate Variants 
#
############################################

#parse geneanno param
refgene="N"
knowngene="N"
ensgene="N"	

if [[ $geneanno =~ "refSeq" ]]
then
	refgene="Y"
fi
if [[ $geneanno =~ "knowngene" ]]
then
	knowngene="Y"
fi
if [[ $geneanno =~ "ensgene" ]]
then
	ensgene="Y"
fi
if [ $hgvs == "N" ]
then
	hgvs=""
fi

#parse verdbsnp/1000g/esp strings
dbsnpstr=${verdbsnp//,/ }
filt_dbsnpstr=${filt_dbsnp//,/ }
g1000str=${ver1000g//,/ }
espstr=${esp//,/ }

if [ $test == "Y" ]
then
	echo "annotate dbsnp: $dbsnpstr"
	echo "annotate esp:   $espstr"
	echo "filter dbsnp: $filt_dbsnpstr"
fi

mutationtaster="N"
avsift="N"
lrt="N"
polyphen2="N"
phylop="N"
ljbsift="N"

#parse old impactscores param (obsolete)
if [[ $impactscores =~ "mutationtaster" ]]
then
	mutationtaster="Y"
fi
if [[ $impactscores =~ "sift" ]]
then
	avsift="Y"
fi
if [[ $impactscores =~ "lrt" ]]
then
	lrt="Y"
fi
if [[ $impactscores =~ "ljbsift" ]]
then
	ljbsift="Y"
fi
if [[ $impactscores =~ "ljb2sift" ]]
then
	ljb2sift="Y"
fi
if [[ $impactscores =~ "pp2" ]]
then
	polyphen2="Y"
fi
if [[ $impactscores =~ "phylop" ]]
then
	phylop="Y"
fi

if [[ $varfile == "Y" ]]
then
	convertcoords="Y"
fi

#ljb refers to Liu, Jian, Boerwinkle paper in Human Mutation with pubmed ID 21520341. Cite this paper if you use the scores

ljb2_sift="N"
ljb2_pp2hdiv="N"
ljb2_pp2hvar="N"
ljb2_lrt="N"
ljb2_mt="N"
ljb2_ma="N"
ljb2_fathmm="N"
ljb2_gerp="N"
ljb2_phylop="N"
ljb2_siphy="N"

# parse ljb2 newimpactscores param
# ljb2_sift, ljb2_pp2hdiv, ljb2_pp2hvar, ljb2_lrt, ljb2_mt, ljb2_ma, ljb2_fathmm, ljb2_gerp++, ljb2_phylop, ljb2_siphy
if [[ $newimpactscores =~ "ljb2_sift" ]]
then
	ljb2_sift="Y"
fi
if [[ $newimpactscores =~ "ljb2_pp2hdiv" ]]
then
	ljb2_pp2hdiv="Y"
fi
if [[ $newimpactscores =~ "ljb2_pp2hvar" ]]
then
	ljb2_pp2hvar="Y"
fi
if [[ $newimpactscores =~ "ljb2_lrt" ]]
then
	ljb2_lrt="Y"
fi
if [[ $newimpactscores =~ "ljb2_mt" ]]
then
	ljb2_mt="Y"
fi
if [[ $newimpactscores =~ "ljb2_ma" ]]
then
	ljb2_ma="Y"
fi
if [[ $newimpactscores =~ "ljb2_fathmm" ]]
then
	ljb2_fathmm="Y"
fi
if [[ $newimpactscores =~ "ljb2_gerp" ]]
then
	ljb2_gerp="Y"
fi
if [[ $newimpactscores =~ "ljb2_phylop" ]]
then
	ljb2_phylop="Y"
fi
if [[ $newimpactscores =~ "ljb2_siphy" ]]
then
	ljb2_siphy="Y"
fi

if [ $otherinfo == "N" ]
then
	otherinfo=""
fi


#column header names we will be adding
# ESP 6500
esp6500si_colheader_ALL="ESP6500si_ALL"
esp6500si_colheader_EA="ESP6500si_EA"
esp6500si_colheader_AA="ESP6500si_AA"
esp6500_colheader_ALL="ESP6500_ALL"
esp6500_colheader_EA="ESP6500_EA"
esp6500_colheader_AA="ESP6500_AA"
esp5400si_colheader_ALL="ESP5400si_ALL"
esp5400si_colheader_EA="ESP5400si_EA"
esp5400si_colheader_AA="ESP5400si_AA"
esp5400_colheader_ALL="ESP5400_ALL"
esp5400_colheader_EA="ESP5400_EA"
esp5400_colheader_AA="ESP5400_AA"


# cg46 cg69
cg46_colheader="CG_46_genomes"
cg69_colheader="CG_69_genomes"

cp $infile originalfile
#run annovar or filter only?
if [ $dorunannovar == "Y" ]
then


	####################################
	#
	#       PREPARE INPUT FILE
	#
	####################################
	
	echo "converting input file"
	vcfheader=""
	if [ $vcf == "Y" ]     #if CG varfile, convert
	then 
		# convert vcf to annovarinput
		$scriptsdir/convert2annovar.pl --format vcf4 --includeinfo --outfile annovarinput $infile 2>&1
		
		#construct header line from vcf file		
		cat $infile | grep "#CHROM" > additionalcols
		sed -i 's/#//g' additionalcols 			
		vcfheader="\t`cat additionalcols`"
		echo "vcfheader:$vcfheader"
		echo -e "chromosome\tbegin\tend\treference\tobserved\t`cat additionalcols`" > originalfile
		cat annovarinput >> originalfile
		
		chrcol=1
		startcol=2
		endcol=3
		refcol=4
		obscol=5

	
	elif [ $varfile == "Y" ]     #if CG varfile, convert
	then 
		# convert varfile
		$scriptsdir/convert2annovar.pl --format cg --outfile annovarinput $infile 2>&1
		echo -e "chromosome\tbegin\tend\treference\talleleSeq\tvarType\thaplotype" > originalfile
		cat annovarinput | cut -f1-6,8 >> originalfile
		cat annovarinput | cut -f1-5 >> annovarinput2
		mv annovarinput2 annovarinput

		chrcol=1
		startcol=2
		endcol=3
		refcol=4
		obscol=5

	elif [ $convertcoords == "Y" ]    # if CG-coordinates, convert 
	then
		#echo "rearranging columns and converting coordinates"
		awk 'BEGIN{
				FS="\t";
				OFS="\t";
			}{
				if(FNR>1) { 
					if( $"'"${vartypecol}"'" == "snp" ){ $"'"${startcol}"'" += 1 }; 
					if( $"'"${vartypecol}"'" == "ins" ){ $"'"${refcol}"'" = "-" };
					if( $"'"${vartypecol}"'" == "del" ){ $"'"${startcol}"'" +=1; $"'"${obscol}"'" = "-" };
					if( $"'"${vartypecol}"'" == "sub" ){ $"'"${startcol}"'" += 1 }; 
			
					printf("%s\t%s\t%s\t%s\t%s\n" ,$"'"${chrcol}"'",$"'"${startcol}"'",$"'"${endcol}"'",$"'"${refcol}"'",$"'"${obscol}"'");				
				}
			}	
			END{
			}' $infile > annovarinput

			#remove any "chr" prefixes
			sed -i 's/chr//g' annovarinput

			awk 'BEGIN{
				FS="\t";
				OFS="\t";
			}{
				if(FNR>=1) { 
					if( $"'"${vartypecol}"'" == "snp" ){ $"'"${startcol}"'" += 1 }; 
					if( $"'"${vartypecol}"'" == "ins" ){ $"'"${refcol}"'" = "-" };
					if( $"'"${vartypecol}"'" == "del" ){ $"'"${startcol}"'" +=1; $"'"${obscol}"'" = "-" };
					if( $"'"${vartypecol}"'" == "sub" ){ $"'"${startcol}"'" += 1 }; 
			
					print $0			
				}
			}	
			END{
			}' $infile > originalfile

			#remove any "chr" prefixes
			sed -i 's/chr//g' originalfile
			sed -i 's/omosome/chromosome/g' originalfile


	else #only rearrange columns if already 1-based coordinates	
		echo "rearranging columns "
		awk 'BEGIN{
				FS="\t";
				OFS="\t";
			}{
				if(FNR>1) { 			
					printf("%s\t%s\t%s\t%s\t%s\n",$"'"${chrcol}"'",$"'"${startcol}"'",$"'"${endcol}"'",$"'"${refcol}"'",$"'"${obscol}"'");				
				}
			}	
			END{
			}' $infile > annovarinput

			#remove any "chr" prefixes
			sed -i 's/chr//g' annovarinput
			sed 's/chr//g' $infile > originalfile
			sed -i 's/omosome/chromosome/g' originalfile
	fi

	echo "...finished conversion"
	



	####################################
	#
	#       RUN ANNOVAR COMMANDS
	#
	####################################

	

	######    gene-based annotation   #######

	# RefSeq Gene
	if [ $refgene == "Y" ]
	then
		echo -e "\nrefSeq gene"
		$scriptsdir/annotate_variation.pl --geneanno --buildver $buildver -dbtype gene ${hgvs} annovarinput $humandb 2>&1
		
		annovarout=annovarinput.variant_function
		sed -i '1i\RefSeq_Func\tRefSeq_Gene\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout  3 4 5 6 7 B.RefSeq_Func,B.RefSeq_Gene

		annovarout=annovarinput.exonic_variant_function 
		sed -i '1i\linenum\tRefSeq_ExonicFunc\tRefSeq_AAChange\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout  
		joinresults originalfile $annovarout  4 5 6 7 8 B.RefSeq_ExonicFunc,B.RefSeq_AAChange
	fi


	# UCSC KnownGene
	if [ $knowngene == "Y" ]
	then
		echo -e "\nUCSC known gene"
		$scriptsdir/annotate_variation.pl --geneanno --buildver $buildver -dbtype knowngene annovarinput $humandb 2>&1
	
		annovarout=annovarinput.variant_function
		sed -i '1i\UCSCKnownGene_Func\tUCSCKnownGene_Gene\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.UCSCKnownGene_Func,B.UCSCKnownGene_Gene

		annovarout=annovarinput.exonic_variant_function
		sed -i '1i\linenum\tUCSCKnownGene_ExonicFunc\tUCSCKnownGene_AAChange\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 4 5 6 7 8 B.UCSCKnownGene_ExonicFunc,B.UCSCKnownGene_AAChange
	fi


	# Emsembl Gene
	if [ $ensgene == "Y" ]
	then
		echo -e "\nEnsembl gene"
		$scriptsdir/annotate_variation.pl --geneanno --buildver $buildver -dbtype ensgene annovarinput $humandb 2>&1
	
		annovarout=annovarinput.variant_function
		sed -i '1i\EnsemblGene_Func\tEnsemblGene_Gene\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.EnsemblGene_Func,B.EnsemblGene_Gene

		annovarout=annovarinput.exonic_variant_function
		sed -i '1i\linenum\tEnsemblGene_ExonicFunc\tEnsemblGene_AAChange\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 4 5 6 7 8 B.EnsemblGene_ExonicFunc,B.EnsemblGene_AAChange
	fi



	######    region-based annotation   #######


	# Transcription Factor Binding Sites Annotation
	if [ $mce == "Y" ]
	then
		echo -e "\nMost Conserved Elements"
	
		if [ $buildver == "hg18" ]
		then
			$scriptsdir/annotate_variation.pl --regionanno --buildver $buildver -dbtype mce44way annovarinput $humandb 2>&1
			annovarout=annovarinput.${buildver}_phastConsElements44way
			sed -i '1i\db\tphastConsElements44way\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
			joinresults originalfile $annovarout 3 4 5 6 7 B.phastConsElements44way	

		else #hg19	
			$scriptsdir/annotate_variation.pl --regionanno --buildver $buildver -dbtype mce46way annovarinput $humandb 2>&1
			annovarout=annovarinput.${buildver}_phastConsElements46way	   
			sed -i '1i\db\tphastConsElements46way\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
			joinresults originalfile $annovarout 3 4 5 6 7 B.phastConsElements46way	
		fi
	
	fi



	# Transcription Factor Binding Sites Annotation
	if [ $tfbs == "Y" ]
	then
		echo -e "\nTranscription Factor Binding Site Annotation"
		$scriptsdir/annotate_variation.pl --regionanno --buildver $buildver -dbtype tfbs annovarinput $humandb 2>&1
	
		# arguments: originalfile, resultfile,chrcol,startcol,endcol,refcol,obscol,selectcolumns
		annovarout=annovarinput.${buildver}_tfbsConsSites
		sed -i '1i\db\tTFBS\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.TFBS	
	fi



	# Identify cytogenetic band for genetic variants
	if [ $cytoband == "Y" ]
	then
		echo -e "\nCytogenic band Annotation"
		$scriptsdir/annotate_variation.pl --regionanno --buildver $buildver -dbtype band annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_cytoBand
		sed -i '1i\db\tBand\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.Band	
	fi


	# Identify variants located in segmental duplications
	if [ $segdup == "Y" ]
	then
		echo -e "\nSegmental Duplications Annotation"
		$scriptsdir/annotate_variation.pl --regionanno --buildver $buildver -dbtype segdup annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_genomicSuperDups
		sed -i '1i\db\tSegDup\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.SegDup	
	fi



	# Identify previously reported structural variants in DGV
	if [ $dgv == "Y" ]
	then
		echo -e "\nDGV Annotation"
		$scriptsdir/annotate_variation.pl --regionanno --buildver $buildver -dbtype dgvMerged annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_dgvMerged
		sed -i '1i\db\tDGV\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.DGV	
	fi


	# Identify variants reported in previously published GWAS studies
	if [ $gwas == "Y" ]
	then
		echo -e "\nGWAS Annotation"
		$scriptsdir/annotate_variation.pl --regionanno --buildver $buildver -dbtype gwascatalog annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_gwasCatalog
		sed -i '1i\db\tGWAS\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.GWAS	
	fi


	
	
	######    filter-based annotation   #######

	#dbSNP
	for version in $dbsnpstr
	do
		if [ $version == "None" ] 
		then
			break
		fi
		echo -e "\ndbSNP region Annotation, version: $version"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype ${version} annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_${version}_dropped		
		sed -i '1i\db\tdb'${version}'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.db${version}	
	

	done



	#1000 Genomes

	if [ $ver1000g != "None" ] 
	then

		for version in $g1000str
		do
			#column headers
			g1000_colheader_ALL="${version}_ALL"
			g1000_colheader_AFR="${version}_AFR"
			g1000_colheader_AMR="${version}_AMR"
			g1000_colheader_ASN="${version}_ASN"
			g1000_colheader_EUR="${version}_EUR"
			g1000_colheader_CEU="${version}_CEU"
			g1000_colheader_YRI="${version}_YRI"
			g1000_colheader_JPTCHB="${version}_JPTCHB"
			
			doALL="N"
			doAMR="N"
			doAFR="N"
			doASN="N"
			doEUR="N"
			doCEU="N"
			doYRI="N"
			doJPTCHB="N"


			if [ $version == "1000g2012apr" ]
			then
				fileID="2012_04"
				doALL="Y"
				if [ $buildver == "hg19" ]
				then
					doAMR="Y"
					doAFR="Y"
					doASN="Y"
					doEUR="Y"
				fi		
			elif [[ $version == "1000g2012feb" && $buildver == "hg19" ]]
			then
				fileID="2012_02"				
				doALL="Y"	
			elif [[ $version == "1000g2010nov" && $buildver == "hg19" ]]
			then
				fileID="2010_11"
				doALL="Y"	
			elif [[ $version == "1000g2010jul" && $buildver == "hg18" ]]
			then
				fileID="2010_07"
				doALL="N"
				doCEU="Y"
				doYRI="Y"
				doJPTCHB="Y"
			else
				echo "unrecognized 1000g version, skipping"
			fi

			#ALL
			if [ $doALL == "Y"  ]
				then
				echo -e "\n1000Genomes ALL"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_all" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_ALL.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_ALL'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_ALL	
			fi

			# AFR
			if [ $doAFR == "Y"  ]
			then
				echo -e "\n1000Genomes AFR"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_afr" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_AFR.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_AFR'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_AFR	
			fi

			
			# AMR
			if [ $doAMR == "Y"  ]
			then
				echo -e "\n1000Genomes AMR"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_amr" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_AMR.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_AMR'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_AMR	
			fi

			# ASN
			if [ $doASN == "Y"  ]
			then
				echo -e "\n1000Genomes ASN"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_asn" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_ASN.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_ASN'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_ASN	
			fi

			# EUR
			if [ $doEUR == "Y"  ]
			then
				echo -e "\n1000Genomes EUR"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_eur" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_EUR.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_EUR'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_EUR	
			fi

			# CEU
			if [ $doCEU == "Y"  ]
			then
				echo -e "\n1000Genomes CEU"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_ceu" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_CEU.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_CEU'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_CEU	
			fi

			# YRI
			if [ $doYRI == "Y"  ]
			then
				echo -e "\n1000Genomes YRI"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_yri" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_YRI.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_YRI'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_YRI	
	

			fi

			#JPTCHB
			if [ $doJPTCHB == "Y"  ]
			then
				echo -e "\n1000Genomes JPTCHB"
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype "${version}_jptchb" annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_JPTCHB.sites.${fileID}_dropped
				sed -i '1i\db\t'$g1000_colheader_JPTCHB'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$g1000_colheader_JPTCHB	
			fi

		done
	fi


	
	
	#### IMPACT SCORE ANNOTATIONS


	if [ $ljb2_sift == "Y" ]
	then
		echo -e "\nLJB2 SIFT Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_sift annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_sift_dropped
		sed -i '1i\db\tLJB2_SIFT\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_SIFT	
	fi
	
	if [ $ljb2_pp2hdiv == "Y" ]
	then
		echo -e "\nLJB2 pp2hdiv Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_pp2hdiv annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_pp2hdiv_dropped
		sed -i '1i\db\tLJB2_PolyPhen2_HDIV\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_PolyPhen2_HDIV	
	fi
	
	if [ $ljb2_pp2hvar == "Y" ]
	then
		echo -e "\nLJB2 pp2hvar Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_pp2hvar annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_pp2hvar_dropped
		sed -i '1i\db\tLJB2_PolyPhen2_HVAR\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_PolyPhen2_HVAR	
	fi
	
	if [ $ljb2_lrt == "Y" ]
	then
		echo -e "\nLJB2 LRT Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_lrt annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_lrt_dropped
		sed -i '1i\db\tLJB2_LRT\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_LRT	
	fi
	
	if [ $ljb2_mt == "Y" ]
	then
		echo -e "\nLJB2 mutationtaster Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_mt annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_mt_dropped
		sed -i '1i\db\tLJB2_MutationTaster\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_MutationTaster	
	fi
	
	if [ $ljb2_ma == "Y" ]
	then
		echo -e "\nLJB2 mutationassessor Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_ma annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_ma_dropped
		sed -i '1i\db\tLJB2_MutationAssessor\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_MutationAssessor	
	fi
	
	if [ $ljb2_fathmm == "Y" ]
	then
		echo -e "\nLJB2 FATHMM Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_fathmm annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_fathmm_dropped
		sed -i '1i\db\tLJB2_FATHMM\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_FATHMM	
	fi
	
	if [ $ljb2_gerp == "Y" ]
	then
		echo -e "\nLJB2 GERP++ Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_gerp++ annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_gerp++_dropped
		sed -i '1i\db\tLJB2_GERP++\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_GERP++	
	fi
	
	if [ $ljb2_phylop == "Y" ]
	then
		echo -e "\nLJB2 PhyloP Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_phylop annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_phylop_dropped
		sed -i '1i\db\tLJB2_PhyloP\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_PhyloP	
	fi
	
	if [ $ljb2_siphy == "Y" ]
	then
		echo -e "\nLJB2 SiPhy Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver $otherinfo -dbtype ljb2_siphy annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb2_siphy_dropped
		sed -i '1i\db\tLJB2_SiPhy\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB2_SiPhy	
	fi



	### OLD IMPACT SCORE ANNOTATIONS

	# SIFT
	if [ $avsift == "Y" ]
	then
		echo -e "\nSIFT Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype avsift annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_avsift_dropped
		sed -i '1i\db\tAVSIFT\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.AVSIFT	
	fi

	#ljb refers to Liu, Jian, Boerwinkle paper in Human Mutation with pubmed ID 21520341. Cite this paper if you use the scores
	# SIFT2
	if [ $ljbsift == "Y" ]
	then
		echo -e "\nLJB SIFT Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype ljb_sift annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb_sift_dropped
		sed -i '1i\db\tLJB_SIFT\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LJB_SIFT
	fi


	# PolyPhen2
	if [ $polyphen2 == "Y" ]
	then
		echo -e "\nPolyPhen Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype ljb_pp2 annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb_pp2_dropped
		sed -i '1i\db\tPolyPhen2\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.PolyPhen2
	fi


	# MutationTaster
	if [ $mutationtaster == "Y" ]
	then
		echo -e "\nMutationTaster Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype ljb_mt annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb_mt_dropped
		sed -i '1i\db\tMutationTaster\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.MutationTaster
	fi


	# LRT
	if [ $lrt == "Y" ]
	then
		echo -e "\nLRT Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype ljb_lrt annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb_lrt_dropped
		sed -i '1i\db\tLikelihoodRatioTestScore\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.LikelihoodRatioTestScore
	fi

	# PhyloP
	if [ $phylop == "Y" ]
	then
		echo -e "\nPhyloP Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype ljb_phylop annovarinput $humandb 2>&1
	
		annovarout=annovarinput.${buildver}_ljb_phylop_dropped
		sed -i '1i\db\tPhyloP\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.PhyloP
	fi


	### ESP  Exome Variant Server
	if [ $esp != "None" ] 
	then
		echo -e "\nESP Annotation"
		for version in $espstr
		do
			echo "version: $version"
			# 6500si ALL
			if [ $version == "esp6500si_all" ]
			then 				
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp6500si_all annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_esp6500si_all_dropped
				sed -i '1i\db\t'$esp6500si_colheader_ALL'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp6500si_colheader_ALL
			fi

			
			# 6500si European American
			if [ $version == "esp6500si_ea" ]
			then
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp6500si_ea annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_esp6500si_ea_dropped
				sed -i '1i\db\t'$esp6500si_colheader_EA'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'" ' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp6500si_colheader_EA
			fi

			# 6500si African Americans
			if [ $version == "esp6500si_aa" ]
			then
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp6500si_aa annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_esp6500si_aa_dropped
				sed -i '1i\db\t'$esp6500si_colheader_AA'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'" ' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp6500si_colheader_AA
			fi


			# 6500 ALL
			if [ $version == "esp6500_all" ]
			then 				
				ls
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp6500_all annovarinput $humandb 2>&1
				
				annovarout=annovarinput.${buildver}_esp6500_all_dropped				
				sed -i '1i\db\t'$esp6500_colheader_ALL'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'" ' $annovarout
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp6500_colheader_ALL
			fi

			
			# 6500 European American
			if [ $version == "esp6500_ea" ]
			then
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp6500_ea annovarinput $humandb 2>&1	
				annovarout=annovarinput.${buildver}_esp6500_ea_dropped
				sed -i '1i\db\t'$esp6500_colheader_EA'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp6500_colheader_EA
			fi

			# 6500 African Americans
			if [ $version == "esp6500_aa" ]
			then
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp6500_aa annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_esp6500_aa_dropped
				sed -i '1i\db\t'$esp6500_colheader_AA'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp6500_colheader_AA
			fi


			# 5400 ALL
			if [ $version == "esp5400_all" ]
			then 				
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp5400_all annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_esp5400_all_dropped
				sed -i '1i\db\t'$esp5400_colheader_ALL'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp5400_colheader_ALL
			fi

			
			# 5400 European American
			if [ $version == "esp5400_ea" ]
			then
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp5400_ea annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_esp5400_ea_dropped
				sed -i '1i\db\t'$esp5400_colheader_EA'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp5400_colheader_EA
			fi

			# 5400 African Americans
			if [ $version == "esp5400_aa" ]
			then
				$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype esp5400_aa annovarinput $humandb 2>&1
	
				annovarout=annovarinput.${buildver}_esp5400_aa_dropped
				sed -i '1i\db\t'$esp5400_colheader_AA'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
				joinresults originalfile $annovarout 3 4 5 6 7 B.$esp5400_colheader_AA
			fi

		done
	fi	


	#GERP++
	if [ $gerp == "Y" ]
	then
		echo -e "\nGERP++ Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype gerp++gt2 annovarinput $humandb 2>&1
	
		annovarout="annovarinput.${buildver}_gerp++gt2_dropped"
		sed -i '1i\db\tGERP++\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.GERP++
	fi


	#COSMIC
	if [[ $cosmic65 == "Y" && $buildver == "hg19" ]]
	then
		echo -e "\nCOSMIC65 Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype cosmic65 annovarinput $humandb 2>&1
	
		annovarout="annovarinput.${buildver}_cosmic65_dropped"
		sed -i '1i\db\tCOSMIC65\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.COSMIC65

	fi

	if [[ $cosmic67 == "Y" && $buildver == "hg19" ]]
	then
		echo -e "\nCOSMIC67 Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype cosmic67 annovarinput $humandb 2>&1
	
		annovarout="annovarinput.${buildver}_cosmic67_dropped"
		sed -i '1i\db\tCOSMIC67\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.COSMIC67

	fi

        if [[ $cosmic70 == "Y" && $buildver == "hg19" ]]
        then
                echo -e "\nCOSMIC70 Annotation"
                $scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype cosmic70 annovarinput $humandb 2>&1

                annovarout="annovarinput.${buildver}_cosmic70_dropped"
                sed -i '1i\db\tCOSMIC70\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout
                joinresults originalfile $annovarout 3 4 5 6 7 B.COSMIC70

        fi	
	if [[ $clinvar == "Y" && $buildver == "hg19" ]]
	then
		echo -e "\nCLINVAR Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype clinvar_20131105 annovarinput $humandb 2>&1
	
		annovarout="annovarinput.${buildver}_clinvar_20131105_dropped"
		sed -i '1i\db\tCLINVAR\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.CLINVAR

	fi
	
	if [[ $nci60 == "Y" && $buildver == "hg19" ]]
	then
		echo -e "\nNCI60 Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype nci60 annovarinput $humandb 2>&1
	
		annovarout="annovarinput.${buildver}_nci60_dropped"
		sed -i '1i\db\tNCI60\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.NCI60

	fi
	
	#cg46
	if [[ $cg46 == "Y"  ]]
	then
		echo -e "\nCG 46 genomes Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype cg46 annovarinput $humandb 2>&1
	
		annovarout="annovarinput.${buildver}_cg46_dropped"
		sed -i '1i\db\t'${cg46_colheader}'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.${cg46_colheader}

	fi


	#cg69
	if [[ $cg69 == "Y"  ]]
	then
		echo -e "\nCG 69 genomes Annotation"
		$scriptsdir/annotate_variation.pl --filter --buildver $buildver -dbtype cg69 annovarinput $humandb 2>&1
	
		annovarout="annovarinput.${buildver}_cg69_dropped"
		sed -i '1i\db\t'${cg69_colheader}'\tchromosome\tstart\tend\treference\talleleSeq"'"$vcfheader"'"' $annovarout 
		joinresults originalfile $annovarout 3 4 5 6 7 B.${cg69_colheader}

	fi


	
	if [ $convertcoords == "Y" ]
	then
		echo "converting back coordinates"
		awk 'BEGIN{
				FS="\t";
				OFS="\t";
			}{
				if (FNR==1)
					print $0
				if(FNR>1) { 
					$"'"${chrcol}"'" = "chr"$"'"${chrcol}"'"
					if( $"'"${vartypecol}"'" == "snp" ){ $"'"${startcol}"'" -= 1 }; 	
					if( $"'"${vartypecol}"'" == "ins" ){ $"'"${refcol}"'" = "" };			
					if( $"'"${vartypecol}"'" == "del" ){ $"'"${startcol}"'" -=1; $"'"${obscol}"'" = "" };
					if( $"'"${vartypecol}"'" == "sub" ){ $"'"${startcol}"'" -= 1 }; 
					print $0
								
				}
			}	
			END{
			}' originalfile > originalfile_coords
	else
		mv originalfile originalfile_coords
	fi

	#restore "chr" prefix?

	#move to outputfile
	if [ ! -s annovarinput.invalid_input ]
	then
		echo "Congrats, your input file contained no invalid lines!" > annovarinput.invalid_input
	fi
	
	cp originalfile_coords $outfile_all
	cp annovarinput.invalid_input $outfile_invalid 2>&1
fi #if $dorunannovar























