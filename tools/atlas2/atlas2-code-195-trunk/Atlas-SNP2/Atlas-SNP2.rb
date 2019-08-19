#!/usr/bin/env ruby

require 'getoptlong'
require 'bigdecimal'

$VERSION="1.4.3"
$REVISION=195

class AtlasSNP2

  #Logistic regression coefficients
    $logistic_FLX={:intercept=>-5.36, :quality_score=>0.11, :NQS_pass=>1.07, :swap=>-4.88, :relpos=>0.91}
    $logistic_XLR={:intercept=>-3.31, :quality_score=>0.1, :NQS_pass=>0.26, :swap=>-3.5, :relpos=>-0.37}
    $logistic_SLX={:intercept=>-9.088, :quality_score=>0.162, :NQS_pass=>1.645, :swap=>0.0, :relpos=>2.349}

  #Sj priori
    $errorPrior_454_less_3=[0.867, 0.867, 0.062, 0.062, 0.063, 0.063, 0.008, 0.008, 0.0003, 0.0003]
    $snpPrior_454_less_3=[0.78, 0.78, 0.122, 0.122, 0.049, 0.049, 0.049, 0.049, 0.00001, 0.00001]

    $errorPrior_454_greater_3=[0.676, 0.676, 0.102, 0.102, 0.074, 0.074, 0.028, 0.028, 0.12, 0.12]
    $snpPrior_454_greater_3=[0.00001, 0.00001, 0.006, 0.006, 0.006, 0.006, 0.041, 0.041, 0.947, 0.947]

    $errorPrior_SLX=[0.991, 0.004, 0.002, 0.001, 0.001, 0.001, 0.000002, 0.000005, 0.0000001, 0.000001]
    $snpPrior_SLX=[0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.027, 0.973]

    $basechange={"A"=>"T","T"=>"A","G"=>"C","C"=>"G","N"=>"N","X"=>"X","*"=>"*","a"=>"t","t"=>"a","g"=>"c","c"=>"g"}
    $swap={0=>"snp",0.5=>"mnp",1=>"swap"}

    SNP = 0
    MNP = 0.5
    SWAP = 1


  def AtlasSNP2.processArguments()

    opts = GetoptLong.new(
      ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
      ["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT],
      ["--help","-h", GetoptLong::NO_ARGUMENT],
      ["--output","-o", GetoptLong::REQUIRED_ARGUMENT],
      ["--target","-t",GetoptLong::OPTIONAL_ARGUMENT],
      ["--maxsub", "-m", GetoptLong::OPTIONAL_ARGUMENT],
      ["--maxindel", "-g", GetoptLong::OPTIONAL_ARGUMENT],
      ["--maxcoverage", "-f", GetoptLong::OPTIONAL_ARGUMENT],
      ["--insertionsize","-p",GetoptLong::OPTIONAL_ARGUMENT],
      ["--always-include","-a",GetoptLong::OPTIONAL_ARGUMENT],
      ["--only-include","-w", GetoptLong::NO_ARGUMENT],
      ["--show-filtered","-F",GetoptLong::OPTIONAL_ARGUMENT],
      ["--sample","-n",GetoptLong::REQUIRED_ARGUMENT],
      ["--cutoff","-c",GetoptLong::OPTIONAL_ARGUMENT],
      ["--mincoverage","-y",GetoptLong::OPTIONAL_ARGUMENT],  
      ["--454_XLR", "-x", GetoptLong::NO_ARGUMENT],
      ["--454_FLX", "-4", GetoptLong::NO_ARGUMENT],
      ["--Illumina", "-s", GetoptLong::NO_ARGUMENT],
      ["--priorerror","-e", GetoptLong::OPTIONAL_ARGUMENT],
      ["--priorerror_454","-l", GetoptLong::OPTIONAL_ARGUMENT]
    )

    optHash = {}
    opts.each do |opt, arg|
      optHash[opt] = arg
    end

    AtlasSNP2.usage() if (optHash.key?("--help"));
    AtlasSNP2.usage() if (optHash.empty?);
    return optHash
  end

  def AtlasSNP2.usage(msg='')
      unless (msg.empty?)
        puts "\n#{msg}\n"
      end
        puts "
Program:  Atlas-SNP2: A SNP callling and evaluation tool for 454 Life Sciences(Roche) FLX/Titanium and Illumina NGS data
Version: #{$VERSION} r#{$REVISION} 
Running Enviroment: Ruby #{RUBY_VERSION}

Usage:	Atlas-SNP2.rb -i [in.sorted.bam] -r [reference.fa] -o [output file] -n [sample name] [choosing platform]

-i	FILE	BAM format alignment file (Required to be sorted by start position)
-r	FILE	FASTA format reference sequence file (Required)
-o	STR	name of output result file (Required)
-n	STR	Sample name used in VCF file (Required)
-t	STR	Only call SNP on given target region (Optional, please refer \"samtools view\" for the target region format)
-a	FILE	file containing sites will always be included(optional)
		-w		only evaluate sites in the list (use with -a, optional)
-F		Include filtered lines in the output that have a QUAL of at least 1

Choosing Platform:
--Illumina 		Illumina
--454_FLX		454 FLX
--454_XLR 		454 Titanium

Setting up prior probability:
-e	FLT	Prior(error|c) when variant coverage number is above 2 for 454 and Illumina data (Default is 0.1)
-l	FLT	Prior(error|c) when variant coverage number is 1 or 2 for 454 data (Default is 0.9)

Setting up filters:
-c      Posterior probability cutoff (Default is 0.95)
-y      Minimal Coverage required for high confidence SNP calls (Default is 6)
-m	FLT	maximum percentage of substitution bases allowed in the alignment (Default is 5.0)
-g	FLT	maximum percentage of insertion and deletion bases allowed in the alignment (Default is 5.0)
-f	INT	maximum number of alignments allowed to be piled up on a site (Default is 1024)
-p	INT	insertion size for pair-end resequencing data (Default is 0)
  "
        exit(2);
  end

  def initialize(optHash)
    @optHash=optHash
    #check ruby version
    if(RUBY_VERSION =~ /(\d+)\.(\d+)/)
    	if($1.to_i < 1 || ($1.to_i == 1 && $2.to_i < 9))
    		STDERR.puts "Atlas-Indel2 requires at least version 1.9 of Ruby, you are running version #{$&}"
    		puts "Quitting"
    		exit 1
    	end
    else
    	raise "failed to detect ruby version"
    end

    #check basic parameters
     if (!@optHash.key?('--input') or !@optHash.key?("--reference") or !@optHash.key?("--output") or !@optHash["--sample"] )
         AtlasSNP2.usage("Option missing!")
         exit(2);
     end

    @SAMFile=@optHash["--input"]
    @outputFile = @optHash["--output"]
    @refFile=@optHash["--reference"]
    @sample = @optHash["--sample"]

    #set parameters

    if @optHash.key?("--target")
	    @target=@optHash["--target"]
    else
	    @target=""
    end

    if @optHash.key?("--maxcoverage")
        @maxcoverage = @optHash["--maxcoverage"].to_i
    else
        @maxcoverage = 1024
    end

    if @optHash.key?("--maxsub")
        @maxsub = @optHash["--maxsub"].to_f
    else
        @maxsub = 5.0
    end

    if @optHash.key?("--maxindel")
        @maxindel = @optHash["--maxindel"].to_f
    else
        @maxindel = 5.0
    end

    if @optHash.key?("--insertionsize")
        @insertionsize = @optHash["--insertionsize"].to_i
    else
        @insertionsize = 0
    end

    #setup priori parameters

    if @optHash.key?("--priorerror")
        @perror_highcov = @optHash["--priorerror"].to_f
        @psnp_highcov = 1 - @perror_highcov
    else
        @perror_highcov = 0.1
        @psnp_highcov = 0.9
    end

    if @optHash.key?("--priorerror_454")
        @perror_lowcov = @optHash["--priorerror_454"].to_f
        @psnp_lowcov = 1 - @perror_lowcov
    else
        @perror_lowcov = 0.9
        @psnp_lowcov = 0.1
    end

    #setup platforms
    
    if @optHash.key?("--454_XLR")
        @platform = "454_Titanium"
        @logistic = $logistic_XLR
        @errPrior0 = $errorPrior_454_greater_3
        @snpPrior0 = $snpPrior_454_greater_3
    elsif @optHash.key?("--Illumina")
        @platform = "Illumina"
        @logistic = $logistic_SLX
        @errPrior0 = $errorPrior_SLX
        @snpPrior0 = $snpPrior_SLX
    elsif @optHash.key?("--454_FLX")
        @platform = "454_FLX"
        @logistic = $logistic_FLX
        @errPrior0 = $errorPrior_454_greater_3
        @snpPrior0 = $snpPrior_454_greater_3
    elsif (!@optHash.key?("--454_FLX") and !@optHash.key?("--Illumina") and !@optHash.key?("--454_XLR"))
	AtlasSNP2.usage("what is the sequencing platform?")
	exit(2)
    end

     puts "#{@platform} data"
     
    if @optHash.key?("--cutoff")
        @cutoff = @optHash["--cutoff"].to_f
    else
        @cutoff = 0.95
    end

    if @optHash.key?("--mincoverage")
        @mincoverage = @optHash["--mincoverage"].to_i
    else
        @mincoverage = 6
    end

    if @optHash.key?("--always-include")
	@vip_sites = Hash.new([])
	File.open(@optHash["--always-include"],'r').each do |line|
		next if line =~ /^#/
		cols = line.split()
		ref, pos = cols[0].sub("chr",""), cols[1].to_i
		@vip_sites[ref] += [pos]
	end
	
	@vip_sites.keys.each do |r|
		@vip_sites[r].sort!
	end
	puts "user specified site list #{@optHash["--always-include"]} loaded"
    end
	
    if @optHash.key?("--show-filtered")
	@show_filtered = true
    else
        @show_filtered = false
    end
   
  end

  def hashReference()
    ref = ''
    @seq = {}
    @snp = {}
    @numdel = {}
    @numref = {}
    refReader=File.open(@refFile,"r")
    refReader.each do |line|
      if line =~ /^>(\S+)/
        ref = $1.sub("chr","")
	puts "loading reference chr#{ref} from \"#{line.chomp}\""
        @seq[ref] = ''
        @numdel[ref] = Hash.new(0)
        @numref[ref] = Hash.new(0)
      else
        @seq[ref] << line.chomp
      end
    end
    refReader.close
   
  end

  def generate_outputfile()
    @outputWriter=File.new("#{@outputFile}.snp","w")
    @outputWriter.print "refName\tcoordinate\trefBase\tvariantBase\toriQual\tvariantReadCov_afterFilter\tnumAlternativeReads_afterFilter\tnumRefReads_afterFilter\ttotalCoverage_afterFilter\tPr(error)j\tPr(SNP)j\tPr(Sj|error,c)\tPr(Sj|SNP,c)\tPrior(error|c)\tPrior(SNP|c)\tPosterior(SNP|Sj,c)j\trefEnv\thomopolymer\treadsInfo\n"
    
    @outputVCF=File.new("#{@outputFile}.vcf","w")
    @outputVCF.puts '##fileformat=VCFv4.0'
    @outputVCF.puts "##fileDate=#{Time.now.strftime("%Y%m%d")}"
    @outputVCF.puts "##source=Atlas-SNP2 v#{$VERSION} r#{$REVISION}"
    @outputVCF.puts "##command=\"#{$command}\""
    @outputVCF.puts "##INFO=<ID=ReqIncl,Number=.,Type=String,Description=\"Site was required to be included in the VCF\">"
    @outputVCF.puts "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    @outputVCF.puts "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
    @outputVCF.puts "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
    @outputVCF.puts "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">"
    @outputVCF.puts "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major Variant Read Depth\">"
    @outputVCF.puts "##FILTER=<ID=low_snpqual,Description=\"SNP posterior probability lower than #{@cutoff}\">"
    @outputVCF.puts "##FILTER=<ID=low_VariantReads,Description=\"Number of variant reads is less than 3\">"
    @outputVCF.puts "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than 0.1\">"
    @outputVCF.puts "##FILTER=<ID=single_strand,Description=\"All variant reads are in a single strand direction\">"
    @outputVCF.puts "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than #{@mincoverage.to_i}\">"
    @outputVCF.puts "##FILTER=<ID=high_coverage,Description=\"Total coverage is more than #{@maxcoverage.to_i}\">"
    @outputVCF.puts "##FILTER=<ID=No_data,Description=\"No valid reads on this site\">"
    @outputVCF.puts "##FILTER=<ID=No_var,Description=\"No valid variant read on this site\">"
    @outputVCF.puts "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t#{@sample}"

  end
  
  def parsing_sam()
      num_filtered = 0
      num_insertionsize = 0
      num_processed = 0
      num_norefence = 0
      num_nm = 0

    IO.popen("samtools view -F 1796 -q 1  #{@SAMFile} #{@target}" ).each do |line|


    begin
	read_sam_cols = line.split("\t")
	query = read_sam_cols[0]
	flag = read_sam_cols[1]
	ref = read_sam_cols[2].sub("chr","")
	offset = read_sam_cols[3].to_i
	mapq = read_sam_cols[4].to_i
	cigar = read_sam_cols[5]
	insertion = read_sam_cols[8].to_i
	query_seq = read_sam_cols[9]
	query_qual = read_sam_cols[10]

	if query_seq.length != query_qual.length
	  $stderr.puts "a broken line"
	  $stderr.puts line
	  next
	end
	
	#apply insertion size as mapping quality control for pair-end data
	if @insertionsize > 0 and (insertion == 0 or insertion > @insertionsize * 4 or insertion < @insertionsize * (-4)) and flag.to_i.to_s(2).reverse[0,1] == '1'
	  num_insertionsize += 1
	  next
	end

	if cigar.count("MIDS") == 0
	  num_unmapped += 1
	  next
	end

	if not @seq.keys.include?(ref)
	  num_norefence += 1
	  next
	end

	#parse optional tags of SAM

	if line =~ /^AS:i:(\S+)$/
	  score = $1
	else
	  score = 'NA'
	end

	if line =~ /NM:i:(\d+)/
	  nm = $1.to_i
	elsif
	  num_nm += 1
	  next
	end

	query_length = query_seq.length


	num_snp = 0
	num_ins = 0
	num_del = 0
	num_clip = 0

	cigar.scan(/([\d]+[MISD])/).each do |part|
	  part.to_s.match(/([\d]+)([MIDS])/)
	  chunk = $1.to_i
	  case $2
	  when "I"
	    num_ins += chunk
	  when "D"
	    num_del += chunk
	  when "S"
	    num_clip += chunk
	  end
	end

	gap = num_ins + num_del
	num_sub = nm - gap
	sub_perc = (num_sub.to_f / (query_length - num_clip) * 10000).round/100.0
	gap_perc = (gap.to_f / (query_length - num_clip) * 10000).round/100.0

	if (sub_perc > @maxsub) or (gap_perc > @maxindel)
	  num_filtered += 1
	  next
	end


	#parse bitwise flag to get strand dir, dir == '1' means reversed strand
	if flag.to_i.to_s(2).reverse[4,1] == '1'
	  dir = '-'
	else
	  dir = '+'
	end

	# parse CIGAR code, find substitutions, calculate qplace, distance_to_3', readenv, tplace, num_snp, num_ins, num_del, cigar_query_length and covered span on reference
	tplace = offset
	qplace = 1
	span = 0
	cigar_query_length = 0
	query_snps = Hash.new([])

	cigar.scan(/([\d]+[MISD])/).each do |part|
	  part.to_s.match(/([\d]+)([MIDS])/)
	  chunk = $1.to_i
	  case $2
	  when "M"
	    cigar_query_length += chunk
	    span += chunk
	      chunk.times do
		begin
		refbase = @seq[ref][tplace - 1,1].upcase
		rescue
		refbase = "N"
		end
		querybase = query_seq[qplace - 1,1].upcase
		
		if num_sub == 0 or refbase == querybase
		  @numref[ref][tplace] += 1
		else
		  allele = query_seq[qplace - 1]
		  qual = query_qual[qplace - 1 ].ord - 33

		  #calculate distance to 3`
		  if dir == '-'
		    dist_3 = qplace - 1
		  else
		    dist_3 = query_length - qplace
		  end

		  nqs = Calculate_NQSpass(qplace - 1, dist_3, query_qual, query_length)

		  #format readEnv
		  if qplace >= 7
		    readenv = query_seq[(qplace - 7)..(qplace - 2)].to_s.downcase + querybase + query_seq[qplace..(qplace + 5)].to_s.downcase
		  elsif qplace == 1
		    readenv = querybase + query_seq[qplace..(qplace + 5)].to_s.downcase
		  else
		    readenv = query_seq[0..(qplace - 2)].to_s.downcase + querybase + query_seq[qplace..(qplace + 5)].to_s.downcase
		  end

		  if dir == '-'
		    readenv.reverse!
		    readenv.gsub!(/[atcgATCGNX*]/) {|s| $basechange[s[0,1]]}
		  end
		  query_snps[tplace] = [allele, qual, qplace, dist_3, nqs, readenv]
		  num_snp += 1
		end

		qplace += 1
		tplace += 1
	      end
	  when "I"
	    qplace += chunk
	    cigar_query_length += chunk
	  when "D"
	    chunk.times do
	      @numdel[ref][tplace] += 1
	      tplace += 1
	    end
	    span += chunk
	  when "S"
	    qplace += chunk
	    cigar_query_length += chunk
	  end
	end
	#All substitutions in this query are stored in query_snps[pos][readinfo]

	if num_snp != num_sub
	  num_nm += 1
	  next
	end

	if cigar_query_length != query_length
	  $stderr.puts "Check CIGAR code failed"
	  $stderr.puts line
	  next
	end
	readinfo = [query, dir, query_length, score, sub_perc, gap_perc]

	#using all substitutions in this query to calculate whether a substitution is a SWAP, SNP or MNP
	query_snps.each_key do |pos|
	  refbase = @seq[ref][pos-1,1].upcase
	  curbase = query_snps[pos][0]
	  next if (refbase == 'N' or curbase == 'N')

	  if query_snps.key?(pos + 1) or query_snps.key?(pos + 2) or query_snps.key?(pos - 1) or query_snps.key?(pos - 2)
	    if query_snps.key?(pos + 1) and refbase == query_snps[pos+1][0] and curbase == @seq[ref][pos,1].upcase
	    query_snps[pos] << SWAP
	    elsif query_snps.key?(pos - 1) and refbase == query_snps[pos-1][0] and curbase == @seq[ref][pos-2,1].upcase
	    query_snps[pos] << SWAP
	    elsif query_snps.key?(pos + 2) and refbase == query_snps[pos+2][0] and curbase == @seq[ref][pos+1,1].upcase
	    query_snps[pos] << SWAP
	    elsif query_snps.key?(pos - 2) and refbase == query_snps[pos-2][0] and curbase == @seq[ref][pos-3,1].upcase
	    query_snps[pos] << SWAP
	    elsif (query_snps.key?(pos + 1) and query_snps[pos + 1][0] != 'N') or (query_snps.key?(pos - 1) and query_snps[pos - 1][0] != 'N')
	    query_snps[pos] << MNP
	    else
	      query_snps[pos] << SNP
	    end
	  else
	    query_snps[pos] << SNP
	  end
	  #complete calculations of all candidate SNPs in this query, pass to a global varible @snp[ref][pos]=[[snp1_readinfo],[snp2_readinfo]..]
	  @snp[ref] = {} unless @snp.key?(ref)
	  @snp[ref][pos] = [] unless @snp[ref].key?(pos)
	  @snp[ref][pos] << query_snps[pos] + readinfo
	end
	
      rescue
	$stderr.puts "a broken line"
	$stderr.puts line
	next
      end

      snp_evaluate(ref, offset)
      num_processed += 1

      if num_processed.modulo(10000) == 0
        print '*'
       print memory_usage = `ps -o rss= -p #{Process.pid}`.to_i
      end

    end

    #Cleanup
    snp_evaluate("ZZZ", 100000000000)
    unless(@vip_sites.nil?) # clear out any VIP sites left (they were on a chromosome that had no data)
    	@vip_sites.keys.each do |r|
	        if @vip_sites[r].length > 0
            vip = @vip_sites[r].shift
            vcfprinter(r, vip, @seq[r][vip - 1, 1], ".", 0, 0, 0, ".","ReqIncl")
          end
        end
    end
    @outputWriter.close

    puts
    puts "#{num_filtered} alignments exceed SUB/INDEL percentage threshold"
    puts "#{num_nm} alignments don't have correct NM tag"
    puts "#{num_insertionsize} alignments don't have proper insertion size for pair-end data"
    puts "#{num_norefence} alignments can not find reference"
    puts "#{num_processed} alignments are processed"

  end

  def snp_evaluate(current_ref, last_pos)
    @snp.keys.sort.each do |ref|
      @snp[ref].keys.sort.each do |pos|
        if pos < last_pos or ref != current_ref
          site_info = "."
          #insert codes to print non-variant VIP sites before pos

          while (!@vip_sites.nil? and @vip_sites[ref].size > 0 and @vip_sites[ref][0] <= pos.to_i)
            vip = @vip_sites[ref].shift
            if @snp[ref][vip].nil?
              vcfprinter(ref, vip, @seq[ref][vip - 1, 1], ".", 0, @numref[ref][vip], 0, ".","ReqIncl")
            else
              site_info = "ReqIncl"
            end
          end
	  
	  next if (@optHash.key?("--only-include") and pos != vip)

          readinfo_output = ''
          added_qual = 0
          refbase = @seq[ref][pos.to_i - 1, 1]
          refbasecov = @numref[ref][pos]


          #format refenv
          if pos.to_i < 7
            refenv = @seq[ref][0, pos.to_i + 6]
          else
            refenv = @seq[ref][pos.to_i - 7, 13]
          end

          alternativereads = @snp[ref][pos].length
          totalcoverage = refbasecov + alternativereads + @numdel[ref][pos]
          next if (totalcoverage > @maxcoverage and pos != vip)
          homonum = homocount(refenv)
          base_num = {}
          @snp[ref][pos].each do |col|
            base=col[0]
            if base_num[base] == nil
              base_num[base] = {}
              base_num[base][:num] = 1
            else
              base_num[base][:num] += 1
            end
          end
          array = base_num.keys.sort {|a,b| base_num[b][:num] <=> base_num[a][:num]}
          varbase = array[0]

          var_cov=0
          rstr=''
          pr_error_j=1
          @snp[ref][pos].each do |read_col|
            env  = read_col[5]
            base = read_col[0]
            qual = read_col[1]
            name = read_col[7]
            score= read_col[10]
            sub  = read_col[11]
            gap  = read_col[12]
            nqs  = read_col[4]
            dist = read_col[3]
            len  = read_col[9]
            dir  = read_col[8]
            swap = read_col[6]
            relpos=dist/len.to_f

            #format readinfo of each substitution
            readinfo_output = "#{base}(#{qual})#{name}(#{dist})(#{score}/#{len})#{dir}#{env}(#{sub}/#{gap})(#{nqs}/#{len})#{$swap[swap]}"

            if base == varbase
                    added_qual += qual
                    #Logistic regression p-value
                    predict = @logistic[:intercept] + @logistic[:quality_score]*qual.to_f + @logistic[:swap]*swap.to_f + @logistic[:NQS_pass]*nqs.to_f + @logistic[:relpos] * relpos.to_f
                    pr_snp_i = Math.exp(predict) / (1 + Math.exp(predict))
                    pr_error_i=1-pr_snp_i
                    pr_error_j*=pr_error_i
                    rstr += readinfo_output+"("+((pr_error_i*1000).round/1000.0).to_s+");"
                    var_cov+=1
            end
          end

            pr_snp_j=1-BigDecimal(pr_error_j.to_s)
            bin = (BigDecimal(pr_snp_j.to_s) *10).floor

            if bin==10
              bin=9
            end

            #set up parameters
            snpPrior = @snpPrior0
            errPrior = @errPrior0
            prior_error_c = @perror_highcov
            prior_snp_c = @psnp_highcov

            if (var_cov < 3) and (@platform == "454_FLX" or @platform == "454_Titanium")
              snpPrior = $snpPrior_454_less_3
              errPrior = $errorPrior_454_less_3
              prior_error_c = @perror_lowcov
              prior_snp_c = @psnp_lowcov
            end

            pr_Sj_error_c = errPrior[bin]
            pr_Sj_SNP_c = snpPrior[bin]
            errorPosterior = pr_Sj_error_c * prior_error_c
            snpPosterior = pr_Sj_SNP_c * prior_snp_c
            pr_SNP_Sj_c_j = 1.0 / (1 + errorPosterior/snpPosterior)

            @outputWriter.print "chr#{ref}\t#{pos}\t#{refbase}\t#{varbase}\t#{added_qual}\t#{var_cov}\t#{alternativereads}\t#{refbasecov}\t#{totalcoverage}\t#{(pr_error_j*1000).round/1000.0}\t#{(pr_snp_j*1000).round/1000.0}\t#{pr_Sj_error_c}\t#{pr_Sj_SNP_c}\t#{prior_error_c}\t#{prior_snp_c}\t#{(pr_SNP_Sj_c_j*1000).round/1000.0}\t#{refenv}\t#{homonum}\t#{rstr}\n"

            vcfprinter(ref, pos, refbase, varbase, var_cov, refbasecov, pr_SNP_Sj_c_j, rstr, site_info)
          end
        end
      #print remaining non-variant required sites before jumping to the next chr

      if (ref != current_ref)
        while (!@vip_sites.nil? and @vip_sites[ref].size > 0)
          vip = @vip_sites[ref].shift
          vcfprinter(ref, vip, @seq[ref][vip - 1, 1], ".", 0, @numref[ref][vip], 0, ".","ReqIncl")
        end
      end
	#clear the @snp @numdel @numref after SNP evaluation
      if ref == current_ref
        next_vip = @vip_sites[ref][0] unless @vip_sites.nil?
        next_vip = 10000000000 if next_vip.nil?
        @snp[ref].delete_if {|key, value| key < last_pos && key < next_vip}
        @numdel[ref].delete_if {|key, value| key < last_pos && key < next_vip}
        @numref[ref].delete_if {|key, value| key < last_pos && key < next_vip}
      else
        @snp.delete(ref)
        @numdel.delete(ref)
        @numref.delete(ref)
      end

    end
  end

  def Calculate_NQSpass(pos, dist, qual, length)
    p_status = 0
    if length > 10
      if dist > length - 6
        p_status  = 1
      elsif pos >= 5 and length - pos >= 6 and qual[pos].ord >= 53 and qual[pos-1].ord >= 48 and qual[pos-2].ord >= 48 and qual[pos-3].ord >= 48 and qual[pos-4].ord >= 48 and qual[pos-5].ord >= 48 and qual[pos+1].ord >= 48 and qual[pos+2].ord >= 48 and qual[pos+3].ord >= 48 and qual[pos+4].ord >= 48 and qual[pos+5].ord >= 48
        p_status = 1
      end
    end
    return p_status
  end

  def homocount(str)
    s = str.size
    lasts = ''
    homo = {}
    st = 0
    flag = 0
    0.upto(s-1) do |i|
      cur = str[i,1].upcase
      if cur == lasts
        if flag == 0 # this is the second
          homo[i-1] = 2
          st = i - 1
        else # this is the third or more
          homo[st] += 1
        end
        flag = 1
      else
        flag = 0
      end
      lasts = cur
    end
    if homo.size > 0
      top = homo.keys.sort {|a,b| homo[b] <=> homo[a]}[0]
      xx = homo[top]
    else
      xx = 1
    end
    return xx
  end

	def passes_strandedness_test(reads_info, num_alt_reads,total_cov )
		return true if total_cov.to_i < 16 #with so few total reads, it is quite possible they are all one direction
		plus_dir = reads_info.scan(/\)([+-])[acgtACGTN]/).count(["+"])
		minus_dir = num_alt_reads.to_i - plus_dir
		return false if [plus_dir.to_f/num_alt_reads.to_f, minus_dir.to_f/num_alt_reads.to_f].min < 0.01
		return true
	end

	def vcfprinter (ref, pos, ref_allele, var_allele, var_num, ref_num, pr, reads_info, site_info)
		#skip low SNP score and non-required sites
		return if (pr < 0.109 and site_info == ".")
			
		
		n = ref_num + var_num
		filter = "" 
		snp_qual = "."
		
		if pr < @cutoff.to_f
			filter += "low_snpqual;"
		end
		
		if n > @maxcoverage
			filter += "high_coverage;"
		end
		
		if n < @mincoverage
			filter += "low_coverage;"
		end

		if !(reads_info == ".") and (! passes_strandedness_test(reads_info, var_num, n))
			filter += "single_strand;"
		end

		if var_num < 3 and var_num > 0
			filter += "low_VariantReads;"
		end

		if var_num > 0
			snp_qual = (-10 * Math.log10(1-pr+0.000001)).round
		end

		if n == 0
			genotype = "."
			filter = "No_data"
		else
			if var_num.to_f/n <= 0.1
				genotype = "0/0"
				filter += "low_VariantRatio"
			elsif var_num.to_f/n > 0.1 and var_num.to_f/n < 0.9
				genotype = "0/1"
			else
				genotype = "1/1"
			end
		end
		
		if (n > 0 and var_num == 0)
			filter = "No_var"
		end
		
		if filter == ""
			filter = "PASS" 
		else
			if(@show_filtered || site_info != ".")
				filter = filter.chomp(";")
			else
				return
			end
		end

		@outputVCF.print ref,"\t",pos,"\t",".","\t",ref_allele,"\t",var_allele,"\t",snp_qual,"\t",filter,"\t",site_info,"\t","GT:VR:RR:DP:GQ","\t","#{genotype}:#{var_num}:#{ref_num}:#{n}:.\n"
	end

end



$command = "#{$0} #{ARGV.join( ' ' )}"
optHash=AtlasSNP2.processArguments()
start=AtlasSNP2.new(optHash)
start.hashReference()
start.generate_outputfile()
start.parsing_sam()

exit(0);
