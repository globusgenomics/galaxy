#!/usr/bin/env ruby

require 'getoptlong'
require 'bigdecimal'


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
      ["--VCF","-v",GetoptLong::OPTIONAL_ARGUMENT],
      ["--sample","-n",GetoptLong::OPTIONAL_ARGUMENT],
      ["--cutoff","-c",GetoptLong::OPTIONAL_ARGUMENT],
      ["--mincoverage","-y",GetoptLong::OPTIONAL_ARGUMENT],  
      ["--XLR", "-x", GetoptLong::NO_ARGUMENT],
      ["--Solexa", "-s", GetoptLong::NO_ARGUMENT],
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
        puts '''
Program:  Atlas-SNP2: A SNP callling and evaluation tool for 454 Life Sciences(Roche) FLX/Titanium and Illumina NGS data
Version:  1.3 (2011-08-18)
Running Enviroment: Ruby 1.9.x

Usage:	Atlas-SNP2.rb -i [in.sorted.bam] -r [reference.fa] -o [output file] [choosing platform] [Setting up VCF output]

-i	FILE	BAM format alignment file (Required to be sorted by start position)
-r	FILE	FASTA format reference sequence file (Required)
-o	STR	    name of output result file (Required)
-t	STR 	Only call SNP on given target region (Optional, please refer "samtools view" for the target region format)

Choosing Platform: (Default is 454 FLX)
-s 		Illumina
-x 		454 Titanium

Setting up VCF output
-v      Output genotypes in VCF format
-n      Sample name used in VCF file (Required when choosing VCF output)
-c      Posterior probability cutoff (Default is 0.95)
-y      Minimal Coverage required for high confidence SNP calls (Default is 8)

Setting up prior probability:
-e	FLT	Prior(error|c) when variant coverage number is above 2 for 454 and Illumina data (Default is 0.1)
-l	FLT	Prior(error|c) when variant coverage number is 1 or 2 for 454 data (Default is 0.9)

Setting up filters:
-m	FLT	maximum percentage of substitution bases allowed in the alignment (Default is 5.0)
-g	FLT	maximum percentage of insertion and deletion bases allowed in the alignment (Default is 5.0)
-f	INT	maximum number of alignments allowed to be piled up on a site (Default is 1024)
-p	INT insertion size for pair-end resequencing data (Default is 0)
  '''
        exit(2);
  end

  def initialize(optHash)
    @optHash=optHash
    
    #check ruby version
    if not RUBY_VERSION =~ /1.9/
        AtlasSNP2.usage("Require ruby 1.9.X or above, your ruby version is #{RUBY_VERSION}")
        exit(2);
    end

    #check basic parameters
     if (!@optHash.key?('--input') or !@optHash.key?("--reference") or !@optHash.key?("--output") )
         AtlasSNP2.usage("Option missing!")
         exit(2);
     end

    #set input files

    @SAMFile=@optHash["--input"]
    @refFile=@optHash["--reference"]
    @target=@optHash["--target"]

    #set parameters

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
    
    if @optHash.key?("--XLR")
        @platform = "454_Titanium"
        @logistic = $logistic_XLR
        @errPrior0 = $errorPrior_454_greater_3
        @snpPrior0 = $snpPrior_454_greater_3
    elsif @optHash.key?("--Solexa")
        @platform = "Illumina"
        @logistic = $logistic_SLX
        @errPrior0 = $errorPrior_SLX
        @snpPrior0 = $snpPrior_SLX
    else
        @platform = "454_FLX"
        @logistic = $logistic_FLX
        @errPrior0 = $errorPrior_454_greater_3
        @snpPrior0 = $snpPrior_454_greater_3
    end
    
    #setup output files
    
    @outputFile = @optHash["--output"]
    
    if @optHash.key?("--VCF")
        if @optHash.key?("--sample")
            @outputVCF = @outputFile + ".vcf"
            @sample = @optHash["--sample"]
        else
            AtlasSNP2.usage("VCF options missing!")
            exit(2);
        end
    end
    
    #setup VCF parameters
    
    if @optHash.key?("--cutoff")
        @cutoff = @optHash["--cutoff"]
    else
        @cutoff = 0.95
    end
    
    if @optHash.key?("--mincoverage")
        @mincoverage = @optHash["--mincoverage"]
    else
        @mincoverage = 8
    end
    

    puts "#{@platform} data"
    
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
        @seq[ref] = ''
        @numdel[ref] = Hash.new(0)
        @numref[ref] = Hash.new(0)
      else
        @seq[ref] << line.chomp
      end
    end
    refReader.close
   puts "Load reference file completed"
  end

  def generate_outputfile()
    @outputWriter=File.new(@outputFile,"w")
    @outputWriter.print "refName\tcoordinate\trefBase\tvariantBase\toriQual\tvariantReadCov_afterFilter\tnumAlternativeReads_afterFilter\tnumRefReads_afterFilter\ttotalCoverage_afterFilter\tPr(error)j\tPr(SNP)j\tPr(Sj|error,c)\tPr(Sj|SNP,c)\tPrior(error|c)\tPrior(SNP|c)\tPosterior(SNP|Sj,c)j\trefEnv\thomopolymer\treadsInfo\n"
  end
  
  def parsing_sam()
      num_filtered = 0
      num_insertionsize = 0
      num_processed = 0
      num_norefence = 0
      num_nm = 0
      
    IO.popen("/usr/bin/samtools view -F 1796 -q 1  #{@SAMFile} #{@target}" ).each do |line|


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
        #complete calculations of all candidate SNPs found in this query, pass them to a global varible @snp[ref][pos]=[[snp1_readinfo],[snp2_readinfo]..]
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
#       print memory_usage = `ps -o rss= -p #{Process.pid}`.to_i
      end

    end
    
    #Cleanup
    snp_evaluate("ZZZ", 100000000000)
    @outputWriter.close
    
    puts "\nProcess SAM files completed"
    puts "#{num_filtered} alignments exceed SUB/INDEL percentage threshold"
    puts "#{num_nm} alignments don't have correct NM tag"
    puts "#{num_insertionsize} alignments don't have proper insertion size for pair-end data"
    puts "#{num_norefence} alignments can not find reference"
    puts "#{num_processed} alignments are processed"
    
    if @optHash.key?("--VCF")
      system("ruby #{__FILE__.sub("Atlas-SNP2.rb", "genotyper.rb")} #{@outputFile} #{@sample} #{@cutoff} #{@mincoverage} > #{@outputVCF}")
    end

  end

  def snp_evaluate(current_ref, last_pos)
    @snp.keys.sort.each do |ref|
      @snp[ref].keys.sort.each do |pos|
        if pos < last_pos or ref != current_ref
          readinfo_output = ''
          added_qual = 0
          refbase=@seq[ref][pos.to_i - 1, 1]

          #format refenv
          if pos.to_i < 7
            refenv = @seq[ref][0, pos.to_i + 6]
          else
            refenv = @seq[ref][pos.to_i - 7, 13]
          end

          refbasecov = @numref[ref][pos]
          alternativereads = @snp[ref][pos].length
          totalcoverage = refbasecov + alternativereads + @numdel[ref][pos]
          next if totalcoverage > @maxcoverage
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
        end
      end
      
      if ref == current_ref
        @snp[ref].delete_if {|key, value| key < last_pos}
        @numdel[ref].delete_if {|key, value| key < last_pos}
        @numref[ref].delete_if {|key, value| key < last_pos}
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

#end class definition
end



optHash=AtlasSNP2.processArguments()
start=AtlasSNP2.new(optHash)
start.hashReference()
start.generate_outputfile()
start.parsing_sam()

exit(0);
