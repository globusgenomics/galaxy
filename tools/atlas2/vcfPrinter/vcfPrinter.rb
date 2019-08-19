#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'getoptlong'
require 'tempfile'
require 'vcf_variant.rb'
require 'sample.rb'
require 'fileutils'
require 'makeSiteFile.rb'
require 'genotype.rb'

#This class reads in per sample vcf files generated from Atlas-SNP, Atlas-Indel and  out multi-sample VCF file

class VcfPrinter

  attr_reader :variants, :samples


  def initialize(vcf_files, out_file, *opt_args)
      @samples = Hash.new
      @variants = Hash.new
      if opt_args.length > 0
        @region = opt_args[0]
      else
        @region = nil
      end
      STDERR.puts Time.now.ctime
      if @region
        siteFileObj = MakeSiteFile.new(vcf_files, out_file, $ARGS['pass'], @region)
      else
        siteFileObj = MakeSiteFile.new(vcf_files, out_file, $ARGS['pass'])
      end
      STDERR.puts 'Sites file generated'
      STDERR.puts Time.now.ctime
      STDOUT.flush
      @samples = siteFileObj.samples
      @variants = siteFileObj.variants
      #puts "VARIANTS:\t#{siteFileObj.variants}"
  end

  def print(outfile)
      #for each sample loop through variants at the end of each sample add a column to the outfile
      if $ARGS['noPileup']
          @samples.each do |sampleName, vcfFile|
              if @region
                sample = Sample.new(sampleName, vcfFile, @region)
              else
                sample = Sample.new(sampleName, vcfFile)
              end
              tmpFile = Tempfile.new('vcfjobs')
              STDOUT.puts "Processing #{sample.name}"
              File.open(outfile,"r").each_line do |line|
                  #puts "#{line}"
                  if line =~ /^#/
                      tmpFile.puts(line)
                  else
                      var = VcfLine.new(line)
                      begin
                          ### ALEX: THIS MIGHT BE THE ISSUE, TRACK WHERE THE GENOTYPE INFO IS GETTING STORED        
                          #puts "#{sample.vcf.get(var.uniqPos, "#{var.chr}_#{var.coor}_#{var.refBase}_.").genotypeInfo}}\n"
                          tmpFile.puts("#{line.strip}\t#{sample.vcf.get(var.uniqPos, "#{var.chr}_#{var.coor}_#{var.refBase}_.").genotypeInfo}")
                          
                      rescue
                          #begin
                          #    tmpFile.puts("#{line.strip}\t#{sample.vcf.get(var.pos, "#{var.chr}_#{var.coor}_.").genotypeInfo}")
                          #rescue
                              tmpFile.puts("#{line.strip}\t.")
                          #end
                      end
                  end
              end
          tmpFile.close()
          FileUtils.mv(tmpFile.path, outfile)
          end
          
      else
          pileups = Hash.new
          $ARGS['pileups'].each do |pileup|
              next if fileDoesNotCheckout(pileup, 'pileup.gz')
              pileups[File.basename(pileup, '.pileup.gz')]=pileup
          end

          @samples.each do |sampleName, vcfFile|
              STDERR.puts Time.now.ctime
              STDOUT.flush
              if @region
                sample = Sample.new(sampleName, vcfFile, @region)
                sample.createPileup(pileups[sampleName], @region)
              else
                sample = Sample.new(sampleName, vcfFile)
                sample.createPileup(pileups[sampleName])
              end
              tmpFile = Tempfile.new('vcfjobs')
              STDOUT.puts "Processing #{sample.name}"
              File.open(outfile,"r").each_line do |line|
                if line =~ /^#/
                  tmpFile.puts(line)
                else
                  var = VcfLine.new(line)
                  begin
                    tmpFile.puts("#{line.strip}\t#{sample.vcf.get(var.uniqPos).genotypeInfo}")
                  rescue
                    begin
                      if ifSNP(var.refBase, var.altBase)
                        tmpFile.puts("#{line.strip}\t#{sample.pileup[var.pos].genotypeInfo('SNP')}")
                      else
                        tmpFile.puts("#{line.strip}\t#{sample.pileup[var.pos].genotypeInfo('INDEL')}")
                      end
                    rescue
                      tmpFile.puts("#{line.strip}\t./.:0:0:0")
                    end
                  end
                end
              end
            tmpFile.close()
            FileUtils.mv(tmpFile.path, outfile)
          end
        end
  end

  def collapseVariants(outFile)
    STDERR.puts "CollapseVariants"
    STDERR.puts Time.now.ctime
    tfile = Tempfile.new('vcfjobs')
    handle = File.open(outFile, 'r')
    handle.each do |line|
      if line =~ /^#/
        tfile.puts line
        next
      end
      vcfObj = VcfLine.new(line)
      new_sample_info = Array.new()
      alt_alleles = [vcfObj.altBase]
      vcfObj.sampleInfoCols.each do |info|
        if info == '.'
          new_sample_info.push(info)
          next
        end
        genoObj = Genotype.new("GT:VR:RR:DP:GQ:FT:AA", info)
        alt_allele = genoObj.get('AA')
        oldGT = genoObj.get('GT')
        if vcfObj.altBase == alt_allele or alt_allele == '.'
          new_sample_info.push(info)
        elsif oldGT =~ /((\.)+\/\.?)|((0)?\/(0))/
          alt_alleles.push(alt_allele) if ! alt_alleles.member? alt_allele
          new_sample_info.push(info)
        else
          alt_alleles.push(alt_allele) if ! alt_alleles.member? alt_allele
          newGT = oldGT.gsub(/[1]/, "#{alt_alleles.index(alt_allele) + 1}")
          new_sample_info.push("#{newGT}:#{info.split(':')[1..-1].join(':')}")
        end
      end
      tfile.puts "#{vcfObj.chr}\t#{vcfObj.coor}\t#{vcfObj.id}\t#{vcfObj.refBase}\
\t#{alt_alleles.join(',')}\t#{vcfObj.qual}\t#{vcfObj.filter}\t.\
\t#{vcfObj.format}\t#{new_sample_info.join("\t")}"
    end
    tfile.close()
    STDERR.puts Time.now.ctime
    FileUtils.mv(tfile.path, outFile)
  end
                                        
  # creates a bed file with list of variant positions                         
  def createBedFile(vcfFile)
      tfile = Tempfile.new('vcfjobs')
      STDERR.puts tfile.path
      `cut -f1,2 #{vcfFile} | grep -v '^#' >#{tfile.path}`
      return tfile.path
  end

  def ifSNP(refAllele, altAllele)
      if refAllele.length==altAllele.length
          return true
      else
          return false
      end
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
      end
  end

  private :createBedFile, :insORdel, :fileDoesNotCheckout, :ifSNP


  def self.usage
    STDERR.puts "USAGE:\n\truby vcfPrinter.rb \n\t -i \"*.vcf\" \n\t -o outputfile \n\t -s region (chr:start-end) \n\t -l \"*.pileup\" \n\t -p <to skip over non PASS variants> \n"
      STDERR.puts "NOTE:\n File naming convention:\n\t [sample_name].vcf <.vcf file> \n\t [sample_name].pileup <.pileup file>"
  end
      

  def self.processArgs()
      if ARGV.length < 4
          STDERR.puts usage()
          Process.exit(0)
      else
          opts = GetoptLong.new(
          ["--vcfFiles","-i",GetoptLong::REQUIRED_ARGUMENT],
          ["--out","-o",GetoptLong::REQUIRED_ARGUMENT],
          ["--region","-s",GetoptLong::OPTIONAL_ARGUMENT],
          ["--pass","-p",GetoptLong::NO_ARGUMENT],
          ["--pileups","-l",GetoptLong::REQUIRED_ARGUMENT],
          ["--usage","-h",GetoptLong::NO_ARGUMENT]
          )
          
          $ARGS = {}
          $ARGS['pass']=false
          $ARGS['noPileup']=true
          $ARGS['region']=false
          
          begin
              opts.each do |opt, arg|
                  case opt
                      when "--vcfFiles"
                          $ARGS['vcf_files'] = Dir.glob(arg).sort 
                      when "--out"
                          $ARGS['out_file'] = arg
                      when "--pass"
                          $ARGS['pass'] = true
                      when "--nopileup"
                          $ARGS['noPileup'] = true
                      when "--pileups"
                          $ARGS['pileups'] = Dir.glob(arg).sort
                          $ARGS['noPileup'] = false
                      when "--region"
                          $ARGS['region'] = arg
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
        tmpOutFile = Tempfile.new('vcfjobsFinalOut')
        puts "TMP FILE: #{tmpOutFile.path}"
        if $ARGS['region']
          vcfPrintObj = VcfPrinter.new($ARGS['vcf_files'], 
                                      tmpOutFile.path, $ARGS['region'])
        else
          vcfPrintObj = VcfPrinter.new($ARGS['vcf_files'],
                                       tmpOutFile.path)
        end
        vcfPrintObj.print(tmpOutFile.path)
        vcfPrintObj.collapseVariants(tmpOutFile.path)
        #FileUtils.mv(tmpOutFile.path, $ARGS['out_file'])
        FileUtils.cp(tmpOutFile.path, $ARGS['out_file'])
        FileUtils.chmod(0644, $ARGS['out_file'])
        STDOUT.puts 'Done!'
end
