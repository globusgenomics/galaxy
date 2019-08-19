#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'getoptlong'
require 'tempfile'
require 'vcf_variant.rb'
require 'sample.rb'
require 'fileutils'
require 'makeSiteFile.rb'
require 'genotype.rb'
require 'functions.rb'

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
              STDERR.puts "Processing #{sample.name}"
              File.open(outfile,"r").each_line do |line|
                  if line =~ /^#/
                      tmpFile.puts(line)
                  else
                      var = VcfLine.new(line)
                      begin        
                          tmpFile.puts("#{line.strip}\t#{sample.vcf.get(var.uniqPos, "#{var.chr}_#{var.coor}_#{var.refBase}_.").genotypeInfo}")
                      rescue
                          tmpFile.puts("#{line.strip}\t.")
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
                      tmpFile.puts("#{line.strip}\t.")
                    end
                  end
                end
              end
            tmpFile.close()
            FileUtils.mv(tmpFile.path, outfile)
          end
        end
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
    STDERR.puts "USAGE:\n\truby vcfPrinter.rb \n\t -i \"*.vcf\" \n\t -o outputfile \n\t -s region (chr:start-end) \n\t -l \"*.pileup\" \n\t --indel Use this flag if input contains Indels \n\t -p <to skip over non PASS variants> \n"
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
          ["--indel",GetoptLong::NO_ARGUMENT],
          ["--usage","-h",GetoptLong::NO_ARGUMENT]
          )
          
          $ARGS = {}
          $ARGS['pass']=false
          $ARGS['noPileup']=true
          $ARGS['region']=false
          $ARGS['indel']=false
          
          begin
              opts.each do |opt, arg|
                  case opt
                      when "--vcfFiles"
                          $ARGS['vcf_files'] = Dir.glob(arg) 
                      when "--out"
                          $ARGS['out_file'] = arg
                      when "--pass"
                          $ARGS['pass'] = true
                      when "--nopileup"
                          $ARGS['noPileup'] = true
                      when "--pileups"
                          $ARGS['pileups'] = Dir.glob(arg)
                          $ARGS['noPileup'] = false
                      when "--region"
                          $ARGS['region'] = arg
                      when "--indel"
                          $ARGS['indel'] = true
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
        if $ARGS['region']
          vcfPrintObj = VcfPrinter.new($ARGS['vcf_files'], 
                                      tmpOutFile.path, $ARGS['region'])
        else
          vcfPrintObj = VcfPrinter.new($ARGS['vcf_files'],
                                       tmpOutFile.path)
        end
        vcfPrintObj.print(tmpOutFile.path)
        if $ARGS['indel']
          Functions.indelProcessing(tmpOutFile.path)
        else
          Functions.collapseVariants(tmpOutFile.path)
        end
        FileUtils.mv(tmpOutFile.path, $ARGS['out_file'])
        FileUtils.chmod(0644, $ARGS['out_file'])
        STDOUT.puts 'Done!'
end
