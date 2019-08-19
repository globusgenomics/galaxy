$:.unshift File.join(File.dirname(__FILE__),'.')

require 'vcf_variant.rb'
require 'tempfile'
require 'fileutils'
require 'genotype'

module Functions
    #Functions involved in pre/post processing vcfs

    def self.commandLineSort(outFile)
        #Using combination if sed,grep and sort to sort the vcf file
        tfile = Tempfile.new('vcfPrinter')
        #substitute X,Y,M with 23,24,25 using sed
        `sed -i -e 's/^X/23/' -e 's/^Y/24/' -e 's/^M/25/' #{outFile}`
        #copy header
        `grep -e '^#' #{outFile} >#{tfile.path}`
        #sort rows
	STDERR.puts 'Sorting sites file...'
        `grep -v -e '^#' #{outFile} | sort -t  $'\t' -k1 -k2 -n >>#{tfile.path}`
        #revert back to X,Y,M from 23,24,25
        `sed -i -e 's/^23/X/' -e 's/^24/Y/' -e 's/^25/M/' #{tfile.path}`
        FileUtils.mv(tfile.path, outFile)
        STDERR.puts 'Sorting done.'
        return outFile
    end

    def self.indelProcessing(outFile)
        STDERR.puts "IndelCleaning"
        STDERR.puts Time.now.ctime
        tfile = Tempfile.new('vcfjobs')
        File.open(outFile,'r').each_line do |line|
          if line =~ /^#/
            tfile.puts line
            next
          end
          new_sample_info = Array.new
          line.strip.split("\t")[9..-1].each do |info|
            new_sample_info.push(info.gsub(/:[A-Z\.]+$/, ''))
          end
          vcfObj = VcfLine.new(line)
          tfile.puts "#{vcfObj.chr}\t#{vcfObj.coor}\t#{vcfObj.id}\t#{vcfObj.refBase}\
\t#{vcfObj.altBase}\t#{vcfObj.qual}\t#{vcfObj.filter}\t.\
\t#{vcfObj.format.chomp(":AA")}\t#{new_sample_info.join("\t")}"
        end
        tfile.close()
        STDERR.puts Time.now.ctime
        FileUtils.mv(tfile.path, outFile)
    end

    def self.collapseVariants(outFile)
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
end
