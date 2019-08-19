#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'getoptlong'
require 'tempfile'
require 'vcf_variant.rb'
require 'sample.rb'
require 'fileutils'
require 'makeSiteFile.rb'

#This class reads in per sample vcf files generated from Atlas-SNP, Atlas-Indel and  out multi-sample VCF file

class VcfPrinter
  def initialize(vcf_files, out_file)
    @samples = Hash.new
    @variants = Hash.new
    puts Time.now.ctime
    siteFileObj = MakeSiteFile.new(vcf_files, out_file, $ARGS['pass'])
    puts Time.now.ctime
    @samples = siteFileObj.samples
    @variants = siteFileObj.variants
    #puts "Number of Samples: #{@samples.length}"
    #puts "Number of Variants: #{@variants.length}"
  end

  def print(outfile)
    #for each sample loop through variants at the end of each sample add a column to the outfile
    if $ARGS['noPileup']
      @samples.each do |sampleName, vcfFile|
        sample = Sample.new(sampleName, vcfFile)
        tmpFile = Tempfile.new('tmp')
        puts "Processing #{sample.name}"
        File.open(outfile,"r").each_line do |line|
          if line =~ /^#/
            tmpFile.puts(line)
          else
            var = VcfLine.new(line)
            begin	
              tmpFile.puts("#{line.strip}\t#{sample.vcf.vcf[var.uniqPos].genotypeInfo}")
            rescue
              tmpFile.puts("#{line.strip}\t./.:.:.:.:.:.")
            end
          end
        end
        tmpFile.close()
        FileUtils.mv(tmpFile.path, outfile)
      end
      
    else
      if $ARGS['hasPileup']
        pileups = Hash.new
        $ARGS['pileups'].each do |pileup|
          next if fileDoesNotCheckout(pileup, 'pileup')
          pileups[File.basename(pileup, '.pileup')]=pileup
        end
        @samples.each do |sampleName, vcfFile|
          sample = Sample.new(sampleName, vcfFile)
          sample.createPileup(pileups[sampleName])
          tmpFile = Tempfile.new('tmp')
          puts "Processing #{sample.name}"
          File.open(outfile,"r").each_line do |line|
            if line =~ /^#/
              tmpFile.puts(line)
            else
              var = VcfLine.new(line)
              #if !sample.vcf.vcf[var.uniqPos].nil?
              begin
                tmpFile.puts("#{line.strip}\t#{sample.vcf.vcf[var.uniqPos].genotypeInfo}")
              rescue
                begin
                  if ifSNP(var.refBase, var.altBase)
                    tmpFile.puts("#{line.strip}\t#{sample.pileup.pileup[var.pos].genotypeInfo('SNP')}")
                  else
                    tmpFile.puts("#{line.strip}\t#{sample.pileup.pileup[var.pos].genotypeInfo('INDEL')}")
                  end
                rescue
                  tmpFile.puts("#{line.strip}\t./.:0:0:0:.:.")
                end
              end
            end
          end
          tmpFile.close()
          FileUtils.mv(tmpFile.path, outfile)
        end
      else
        bams = Hash.new
        $ARGS['bams'].each do |bam|
          next if fileDoesNotCheckout(bam, 'bam')
          bams[File.basename(bam, '.bam')]=bam
        end
        bedFile = createBedFile(outfile)
        @samples.each do |sampleName, vcfFile|
          bamFile = bams[sampleName]
          sample = Sample.new(sampleName, vcfFile)
          sample.generatePileup(bamFile, bedFile, $ARGS['reference'])
          tmpFile = Tempfile.new('tmp')
          puts "Processing #{sample.name}"
          File.open(outfile,"r").each_line do |line|
            if line =~ /^#/
              tmpFile.puts(line)
            else
              var = VcfLine.new(line)
              #tmpFile.puts("#{line.strip}\t#{sample.vcf.vcf[var.uniqPos].genotypeInfo}") if !sample.vcf.vcf[var.uniqPos].nil?
              begin
                tmpFile.puts("#{line.strip}\t#{sample.vcf.vcf[var.uniqPos].genotypeInfo}")
              rescue
                begin
                  if ifSNP(var.refBase, var.altBase)
                    tmpFile.puts("#{line.strip}\t#{sample.pileup.pileup[var.pos].genotypeInfo('SNP')}")
                  else
                    tmpFile.puts("#{line.strip}\t#{sample.pileup.pileup[var.pos].genotypeInfo('INDEL')}")
                  end
                rescue
                  tmpFile.puts("#{line.strip}\t./.:0:0:0:.:.")
                end
              end
            end
          end
          tmpFile.close()
          FileUtils.mv(tmpFile.path, outfile)
        end
        bedFile.unlink
      end
    end
  end
					

  #creates INFO column in the final multi-sample VCF
  def generateInfoField(outFile)
    tFile = Tempfile.new('tmp')
    puts tFile.path
    File.open(outFile,"r").each_line do |line|
      if line =~ /^#/
        tFile.puts(line)
      else
        totalDepth = 0
        totalSample = 0
        line = line.strip.split("\t")
        (9..line.length-1).each do |i|
          depth = line[i].split(':')[3]
          next if line[i] == './.:.:.:.:.:.' or line[i] == './.:0:0:0:.:.'
          totalSample += 1
          next if depth == '.' or depth == '0'
          totalDepth += depth.to_i
        end
        infoField = "NS=#{totalSample};DP=#{totalDepth}"
        
        #modified by liubo, because of the following error 
        #./vcfPrinter/vcfPrinter.rb:147:in `block in generateInfoField': undefined method `join' for nil:NilClass (NoMethodError)
        
        #tFile.puts("#{line[0..6].join("\t")}\t#{infoField}\t#{line[8..line.length-1].join("\t")}")
        tFile.puts("#{line[0..6]}\t#{infoField}\t#{line[8..line.length-1]}\t")
      end
    end
    tFile.close()
    FileUtils.mv(tFile.path, outFile)
  end
  
  
  # collapses altBase
  def collapseVariants(outFile)
    tfile = Tempfile.new('tmp')
    file = File.open(outFile, 'r')
    prevLine = ''
    while (currentLine=file.gets)
      if currentLine =~ /^#/
        tfile.puts(currentLine)
      else
        if prevLine == ''
          prevLine = currentLine
        else
          prevObj = VcfLine.new(prevLine)
          currentObj = VcfLine.new(currentLine)
          #if the chr and coor matches previous line then they need to collapsed
          if prevObj.chr == currentObj.chr and prevObj.coor == currentObj.coor
            #update the current line with alt and ref allele info from previous line
            currentObj.update(prevObj)
            currentGenotypes=currentObj.sampleInfoCols()
            prevGenotypes=prevObj.sampleInfoCols()
            #make two bins puts genotypes you want to change in one and don't inthe other
            genoOfInterest={change:{}, noChange:{}}
            #genotypes in the current line need to be changed so add them to the change bin
            currentGenotypes.each do |x|
              if x.split(':')[0] =~ /((\d)+\/\.?)|((\.)?\/(\d))|((\d)+\/(\d)+)/
                genoOfInterest[:change].store(x,currentGenotypes.index(x))
              else
                next
              end
            end
            #genotypes from the previous line dont change so add them to noChange bin
            prevGenotypes.each do |y|
              #if the genotype feild looks like 1/1,1/0,0/1,1/.,./1 then add it to the bin
              if y.split(':')[0] =~ /((\d)+\/\.?)|((\.)?\/(\d))|((\d)+\/(\d)+)/
                genoOfInterest[:noChange].store(prevGenotypes.index(y),y)
              else
                next
              end
            end
            #this method will update genotypes in vcfObj
            currentObj.updateSampleInfo(genoOfInterest)
            pline = "#{currentObj.chr}\t#{currentObj.coor}\t#{currentObj.id}\t#{currentObj.refBase}\t#{currentObj.altBase}\t#{currentObj.qual}\t#{currentObj.filter}\t#{currentObj.info}\t#{currentObj.format}\t#{currentObj.sampleInfo.join("\t")}"
            prevLine=pline
          else
            tfile.puts(prevLine)
            prevLine = currentLine
          end
        end
      end
    end
    tfile.puts(prevLine)
    tfile.close()
    FileUtils.mv(tfile.path, outFile)
  end
  
  def collapseComplexCases(outFile)
    tfile = Tempfile.new('tmp')
    file = File.open(outFile, 'r')
    while (currentLine=file.gets)
      if currentLine =~ /^#/
        tfile.puts(currentLine)
      else
        #checking to see if there are multiple ref alleles
        if currentLine.split("\t")[3].split(',').length > 1
          #ref.zip(alt).map(insORdel)=>[ref{tab}alt{tab}INS]
          ref_alt_zipd = currentLine.split("\t")[3].split(',').zip(currentLine.split("\t")[4].split(',')).map {|ref,alt| "#{ref}\t#{alt}\t#{insORdel(ref,alt)}"}
          #get me the largest ref allele
          bigRef=currentLine.split("\t")[3].split(',').inject do |x,y| x.length > y.length ? x:y end
          #container for new modified alleles
          newAltAlleles=[]
          ref_alt_zipd.each do |allele|
            allele=allele.split("\t")
            if allele[2]=='INS'
              #if its an insertion the compare the altallele to big ref allele and stick the ins bases after the first bigRef base.
              newAllele="#{bigRef[0]}#{allele[1][1,allele[1].length]}#{bigRef[1,bigRef.length]}"
              newAltAlleles.push(newAllele)
            elsif allele[2]=='DEL'
              #if its a del then if is the bigReg then add altAllele to container if not then find the del bases in bigRef and substitute
              if allele[0]==bigRef
                newAltAlleles.push(allele[1])
              else
                newAllele=bigRef.sub(allele[0][1,allele[0].length],'')
                newAltAlleles.push(newAllele)
              end
            elsif allele[2]=='SNP'
              if allele[0]==bigRef
                newAltAlleles.push(allele[1])
              else
                puts "SNP error in collapseComplexCases #{allele.join(',')}"
                Process.exit(0)
              end
            else
              next
            end
          end
          currentObj=VcfLine.new(currentLine)
          pline = "#{currentObj.chr}\t#{currentObj.coor}\t#{currentObj.id}\t#{bigRef}\t#{newAltAlleles.join(",")}\t#{currentObj.qual}\t#{currentObj.filter}\t#{currentObj.info}\t#{currentObj.format}\t#{currentObj.sampleInfoCols.join("\t")}"
          tfile.puts(pline)
        else
          tfile.puts(currentLine)
        end
      end
    end
    tfile.close()
    FileUtils.mv(tfile.path, outFile)
  end
  
  
  # creates a bed file with list of variant positions 			
  def createBedFile(vcfFile)
    random_string = randomFileNameSuffix(7)
    filename = "/tmp/#{random_string}"
    File.open(filename, 'w') do |tfile|
      tfile.puts(`cut -f1,2 #{vcfFile} | grep -v '^#'`)
    end
    puts filename
    #tfile = Tempfile.new('tmp')
    #`cut -f1,2 #{vcfFile} | grep -v '^#' >#{tfile.path}`
    #puts "cut -f1,2 #{vcfFile} | grep -v '^#' >#{tfile.path}"
    #return tfile.path
    return filename
  end

  def randomFileNameSuffix (numberOfRandomchars)
    s = ""
    numberOfRandomchars.times { s << (65 + rand(26))  }
    s
  end  

  def insORdel(refAllele, altAllele)
    if refAllele.length > altAllele.length
      return "DEL"
    elsif refAllele.length < altAllele.length
      return "INS"
    else
      return 'SNP'
    end
  end
  
  def fileDoesNotCheckout(file, fileType)
    if File.exist?(file) and File.fnmatch("*.#{fileType}",file)
      return false
    else
      raise "Either #{file} does not exist or the file does not have proper .#{fileType} extension"
      return true
      #return false
    end
  end
  
  
  def ifSNP(refAllele, altAllele)
    if refAllele.length==altAllele.length
      return true
    else
      return false
    end
  end
  
  private :createBedFile, :insORdel, :fileDoesNotCheckout, :ifSNP
  
  def self.usage
    puts "USAGE:\n\truby vcfPrinter.rb \n\t -i \"*.vcf\" \n\t -o outputfile \n\t -b \"*.bam\" \n\t -r reference.fasta \n\t -l \"*.pileup\" \n\t -p <to skip over non PASS variants> \n\t -n <This flag will lead to no pileup generated info in VCF file> \n"
    puts "NOTE:\n File naming convention:\n\t [sample_name].vcf <.vcf file> \n\t [sample_name].bam <.bam file> \n\t [sample_name].pileup <.pileup file>"
  end
  
  
  def self.processArgs()
    if ARGV.length < 4
      puts usage()
      Process.exit(0)
    else
      opts = GetoptLong.new(
                            ["--vcfFiles","-i",GetoptLong::REQUIRED_ARGUMENT],
                            ["--out","-o",GetoptLong::REQUIRED_ARGUMENT],
                            ["--bams","-b",GetoptLong::REQUIRED_ARGUMENT],
                            ["--reference","-r",GetoptLong::REQUIRED_ARGUMENT],
                            ["--pass","-p",GetoptLong::NO_ARGUMENT],
                            ["--pileups","-l",GetoptLong::REQUIRED_ARGUMENT],
                            ["--nopileup","-n",GetoptLong::NO_ARGUMENT],
                            ["--usage","-h",GetoptLong::NO_ARGUMENT]
                            )
      $ARGS = {}
      $ARGS['pass']=false
      $ARGS['noPileup']=false
      $ARGS['hasPileup']=false
      begin
        opts.each do |opt, arg|
          case opt
          when "--vcfFiles"
            $ARGS['vcf_files'] = Dir.glob(arg) 
          when "--out"
            $ARGS['out_file'] = arg
          when "--bams"
            $ARGS['bams'] = Dir.glob(arg)
          when "--reference"
            $ARGS['reference'] = arg
          when "--pass"
            $ARGS['pass'] = true
          when "--nopileup"
            $ARGS['noPileup'] = true
          when "--pileups"
            $ARGS['pileups'] = Dir.glob(arg)
            $ARGS['hasPileup'] = true
          when "--usage"
            puts usage()
            Process.exit(0)
          end
        end
        return $ARGS
      rescue
        VcfPrinter.usage()
        Process.exit(0)
      end
    end	
  end
  
end

if __FILE__ == $0
  VcfPrinter.processArgs()
  #siteFileObj = MakeSiteFile($ARGS['vcf_files'], $ARGS['out_file'], $ARGS['pass']) 
  #	tmpOutFile = Tempfile.new('vcfjobsFinalOut', '/space1/tmp')   #modified by liubo, because there's no /space1 directory
  tmpOutFile = Tempfile.new('vcfjobsFinalOut', '/tmp')
  vcfPrintObj = VcfPrinter.new($ARGS['vcf_files'], tmpOutFile.path)
  vcfPrintObj.print(tmpOutFile.path)
  vcfPrintObj.collapseVariants(tmpOutFile.path)
  vcfPrintObj.generateInfoField(tmpOutFile.path)
  vcfPrintObj.collapseComplexCases(tmpOutFile.path)
  FileUtils.mv(tmpOutFile.path, $ARGS['out_file'])
  FileUtils.chmod(0644, $ARGS['out_file'])
  puts 'Done!'
end
