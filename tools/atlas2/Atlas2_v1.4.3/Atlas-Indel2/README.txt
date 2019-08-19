Atlas-Indel2 for Indel Calling in Whole Exome Capture Sequencing Data
version 1.4.1 Sept 7, 2012
AUTHORS: Danny Challis, Uday Evani, Fuli Yu, Human Genome Sequencing Center, Baylor College of Medicine
CONTACT: challis@bcm.edu, fyu@bcm.edu

Atlas-Indel2 is a suite of indel calling tools based on logistic regression models 
made up of pertinent variables sequencing data.  This version of Atlas-Indel includes
regression models for Illumina and SOLiD exome data.


SYSTEM REQUIREMENTS:

    * Unix-like operation systems
    * Ruby 1.9.1: http://www.ruby-lang.org/en/downloads/
    * SAMtools must be installed and runnable by invoking the "samtools" command
      -SAMtools may be obtained at http://samtools.sourceforge.net/
    

LICENSE:
Copyright (c) 2012, Human Genome Sequencing Center, Baylor College of Medicine

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

	*Redistributions of source code must retain the above copyright notice, this
	 list of conditions and the following disclaimer.

	*Redistributions in binary form must reproduce the above copyright notice, this
	 list of conditions and the following disclaimer in the documentation and/or
	 other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE


DATA PRE/POST PROCESSING REQUIREMENTS: 
We recommend local realignment around indels and high variation regions using GATK or 
other third party tools.  While you may run Atlas-Indel2 without local realignment, greater 
sensitivity is possible with it. Details on local realignment using GATK can be found at 
http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels.

We do not recommend base quality recalibration for SOLiD data.  The recalibration process
significantly lowers base qualities which severly limits sensitivity.

Atlas-Indel2 now includes a simple genotyping algorithm.  Future versions will include a more 
robust and accurate genotyping model.

example pipeline:
	*map fastq files using BWA to get BAM files
	*sort BAM files
	*locally realign BAMs using GATK local realigner
	*run Atlas-Indel2 on BAMs to get individual VCF files
	*merge VCF files using vcfPrinter (included)


USAGE:
	ruby Atlas-Indel2.rb -b input_bam -r reference_sequence -o outfile [-S or -I]
		-S --solid-exome-model (use the model trained for SOLiD exome data)
		-I --illumina-exome-model (use the model trained for Illumina exome data)

		optional arguments:
		-s --sample [by default taken from infile name] (the name of the sample)
		-p --p-cutoff (the p cutoff for the regression model)
		-P --p-1bp-cutoff (stricter p cutoff for 1bp indels)
		-B --bed (bed file to specify regions to be called)
		-a --always-include (file of sites with annotation to always include in VCF)
                -F --show-filtered (include filtered indels with a QUAL>=1 in the VCF)
		-O --orig-base-qual
		-N --norm-base-qual
		-t --min-total-depth
		-m --min-var-reads 
		-v --min-var-ratio
		-f --strand-dir-filter (requires at least one read in each direction, 
					extremely limits sensitivity)
		-n --near-read-end-ratio (filters indels where more than the specified 
					fraction of indel reads are near the read end)
		-h --homo-var-cutoff' (homozygous variant cutoff for genotyping)


This usage information can be viewed at any time by running the program without arguments.

Mandatory arguments:

-b, --bam=FILE
	The input BAM file.  It must be sorted. It does not need to be
	indexed.  A read mask of 1796 is used on the bitwise flag.

-r, --reference=FILE
	The reference sequence to be used in FASTA format.  This must be the same
	version used in mapping the sequence. It does not need to be indexed.

-o, --outfile=FILENAME
	The name of the output VCF file.  Output is a bare-bones VCFv4 file with
	a single sample.  These files can be merged into a more complete VCF
	file using the vcfPrinter (included).  If the file already exists, it will
	be overwritten. NOTE: For use with vcfPrinter, you should name your vcf file
	the same as your BAM file, simply replacing ".bam" with ".vcf"

-I or -S
	You must include one of these flags to specify either the Illumina or SOLiD
	regression model to be used.

Optional arguments:
Note: Different platform models have different defaults.

-s, --sample=STRING 
	The name of the sample to be listed in the output VCF file.  If not specified
	the sample name is harvested from the input BAM file name, taking the
	first group of characters before a . (dot) is found.  For example, with
	the filename "NA12275.chrom1.bam" the sample name would be "NA12275".

-p, --p-cutoff=FLOAT
	Defaults: Illumina:0.5, Solid:0.5
	The indel probability (p) cutoff value for the logistic regression model.  Indels with a 
	p less than this cutoff will not be called.  Increasing this cutoff will
	increase specificity, but will lower sensitivity. If you adjust this cutoff,
	you should usually also adjust the p-1bp-cutoff (see below).

-P --p-1bp-cutoff
	Defaults: Illumina:0.5, Solid:0.88
	The indel probability (p) cutoff value for 1bp indels.  This may be set to a stricter standard
	than the normal p-cutoff to increase callset specificity.  This is very useful 
	for SOLiD data, but should not be generally needed for Illumina data.

-B --bed=FILE
	Here you may specify a bed file which contians the region you wish to limit
	your indel calling to.  Only reads inside the region will be process, 
	which can significantly shorten the runtime.

-a --always-include=FILE
	Here you may specify a tab-delimited file of genome positions which
	will always be included in the output VCF regardless of if there is an
	indel or even data there.  
	
-F --show-filtered
	This option will cause the filtered indels to also be included in the
	output VCF and marked with explanitory filters. In order to keep the output 
	to a reasonable size, only indels with a QUAL<1 are excluded.

-O --orig-base-qual
	This is the default for SOLiD, it is not recommended for Illumina data.
	This option has the algorithm use the original base qualities as specified in
	the OQ tag if included in the BAM file.  If the BAM file does not include OQ
	tags, the normal base quality is used.

-N --norm-base-qual
	This is the default for Illumina, it is not recommended for SOLiD data.
	This option specifies the algorithm should use the normal base qualities,
	as specified in the QUAL column of the BAM file.


Heuristic Cutoffs:
Most of these variables have already been considered by the regression model, so you
shouldn't usually need to alter them.  However you are free to change them to meet 
your specific project requirements.

-t, --min-total-depth=INT
	Defaults: Illumina:5, SOLiD:2
	The minimum total depth coverage required at an indel site.  Indels at a
	site with less depth coverage will not be called. Increasing this value will 
	increase specificity, but lower sensitivity.
	Suggested range: 2-12

-m, --min-var-reads=INT
	Defaults: Illumina:2, SOLiD:2
	The minimum number of variant reads required for an indel to be called.
	Increasing this number may increase specificity but will lower sensitivity.
	Suggested range: 1-5

-v, --min-var-ratio=FLOAT
	Defaults: Illumina:0.06, SOLiD:0.05
	The variant-reads/total-reads cutoff.  Indels with a ratio less than the
	specified value will not be called.  Increasing this value may increase
	specificity, but will lower sensitivity. 
	Suggested range: 0-0.1

-f, --strand-dir-filter
	Default: Illumina:disabled, SOLiD:disabled
	When included, requires indels to have at least one variant read in each
	strand direction.  This filter is effective at increasing the
	specificity, but also carries a heavy sensitivity cost.

-n, --near-read_end_ratio=FLOAT
	Default: Illumin:0.8, SOLiD:1.0 (disabled)
	The read end ratio is defined as the number of variant reads where the 
	variant is within 5bp of a read end divided by the total variant read
	depth.  If this ratio is greater than the specified value, the indel
	is filtered.
	Suggested range: 0.7-1.0
	
-h, --homo-var-cutoff 
	Default: Illumina:0.6, Solid:0.5
	The homozygous variant cutoff.  This cutoff is used in the preliminary 
	genotyping performed by Atlas-Indel2.  If the variant reads divided by
	the variant reads + the reference reads is greater than this cutoff it
	will be marked as a homozygote, otherwise it will be marked as a 
	heterozygote.

EXAMPLES:
ruby Atlas-Indel2.rb -b NA12275.chrom1.bam -r ~/refs/human_g1k_v37.fasta -o ~/NA12275.chrom1.vcf -p 0.55 -I
ruby Atlas-Indel2.rb -b seq1.10.2010.chrom1.bam -r ~/refs/human_g1k_v37.fasta -o ~/NA12275.chrom1.vcf -t 10 -m 5 -s NA12275 -B target_region.bed -S


CHANGES:
v1.4.1
* Added alway-include option
* Added show-filtered option
* Fixed bug caused by passing a non-fasta reference genome
* Fixed bug occationally returning infinite P value in INFO column
* Fixed bug caused by reads mapping past the end of the reference genome
* Made Atlas-Indel2 more tolerant of malformed SAM lines

v1.0
* updated SOLiD model and adjusted P cutoffs
* changed -P cutoff to apply to both 1bp insertions and deltions (rather than
  just 1bp deletions)


v0.3.1
* added options to use original base quality
* fixed bug that sometimes returned success error code after a failure
* fixed bug in simple_genotyper that caused samples with exactly 0.05 variant read ratio to be 0/0
* fixed bug in simple genotyper that caused genotypes to occationaly read ./.
* fixed bug in bed_filter that was filtering some on-target reads in very small target regions

v0.3
* updated SOLiD and Illumina models and recalibrated default settings
* Implemented the ability to input a bed file to call only on-target indels
* switched from using z cutoffs to using p cutoffs
* modified 1bp p cutoff to only filter 1bp deletions
* fixed bug where the strand direction filter failed to be enabled
* Added check for proper ruby version
* fixed bug that occasionally allows an indel quality of 110 (max should be 100)
* minor code-structure changes

v0.2.1
* added read_level model and improved site level model for SOLiD data
* adjusted default SOLiD z cutoff to 0.0 (to reflect new model)
* added check for proper ruby version
* minor codes structure changes
* added additional heuristic filter that allows for a stricter z cutoff for 1bp indels, very useful for SOLiD data
* integrated heuristic genotyping -implemented
* fixed bug where Atlas-Indel2 crashes if a BAM chromosome is not in the reference
* now will keep ‘chr’ in the chromosome label if it is in the BAM
* the depreciated script "Atlas-Indel2-Illum-Exome.rb, has been removed.  Please use Atlas-Indel2.rb with the -I flag instead.

v0.2
* Implemented regression model for SOLiD data.  You must now specify a regression model (-S or -I).
* Renamed main script to Atlas-Indel.rb.
* Modified Reference sequence class to allow for unsorted reference genomes.
* Added the indel z to the info column of the VCF output (not included after running VCF printer).
* Now echos all settings back onto the command line.
* Fixed a bug that caused loss of precision in the normalized variant square variable of the Illumina site model.
* Fixed a bug in the depth coverage algorithm that caused reads not to be counted in total depth at the deleted sites.
* Fixed the sample columns order to be comaptible with vcfPrinter.
* Removed "x flagged lines skipped" message at end of run.
