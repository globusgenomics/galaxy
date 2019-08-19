
$:.unshift File.join(File.dirname(__FILE__),'.')

require 'settings.rb'
require 'vcf.rb'
require 'vcf_variant.rb'
require 'tempfile'

class MakeSiteFile
	attr_reader :samples, :variants

        def initialize(vcf_files, out_file, onlyPASS=false, *opt_args)
          @samples = Hash.new
          @variants = Hash.new
	  abort("Did not find vcf files. Trying giving the complete path to the location of vcf files.\n Terminating process.") if vcf_files.length==0
          vcf_files.each do |vcf_file|
            next if fileDoesNotCheckout(vcf_file)
            @samples[File.basename(vcf_file, '.vcf.gz')]=vcf_file
          end
          STDERR.puts "Number of Samples: #{@samples.length}"

          #Check if region is specified else process complete vcf file
          if opt_args.length > 0
            @region = opt_args[0]
          else
            @region = false
          end

          #Loop through all VCF Files and collect all variants in a large Hash
          @samples.each do |sampleName, fileName|
            if @region
              vcfObj = Vcf.new(fileName, @region)
            else
              vcfObj = Vcf.new(fileName)
            end
            vcfObj.each do |key, var|
	      next if var.altBase=='.'
              if var.altBase=='.'
                  if !@variants[var.pos].nil?
                      @variants[var.pos].updateQual(var.qual)
                      @variants[var.pos].filter = var.filter
                  else
                      @variants[var.pos] = var
                  end
              else
                if onlyPASS
                  next if var.filter!='PASS'
                  if !@variants[var.uniqPos].nil?
                    @variants[var.uniqPos].updateQual(var.qual)
                    @variants[var.uniqPos].filter = var.filter
                  else
                    @variants[var.uniqPos] = var
                  end
                else
                  if !@variants[var.uniqPos].nil?
                    @variants[var.uniqPos].updateQual(var.qual)
                    @variants[var.uniqPos].filter = var.filter
                  else
                    @variants[var.uniqPos] = var
                  end
                end
              end
            end
          end

          STDERR.puts "Number of Sites: #{@variants.length}"
          #STDOUT.puts "Analyzing #{chr}:#{start_pos}-#{end_pos}"

          File.open(out_file,"w") do |writer|
            #Printing VCF header information
            writer.puts "##fileformat=VCFv4.0"
            time = Time.new
            writer.puts "##fileDate=#{time.ctime}"
            writer.puts "##source=VcfPrinter v#{VERSION} r#{REVISION}"
            #writer.puts "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
            #writer.puts "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
            writer.puts "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
            writer.puts "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major variant Read Depth\">"
            writer.puts "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">"
            writer.puts "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">"
            writer.puts "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
            writer.puts '##FORMAT=<ID=FT,Number=.,Type=String,Description="Sample Genotype Filter">'
            writer.puts '##FORMAT=<ID=AA,Number=1,Type=String,Description="Alt allele associated with variant depth">'
            writer.puts "##FILTER=<ID=low_snpqual,Description=\"SNP posterior probability lower than 0.95 (default)\">"
            writer.puts "##FILTER=<ID=low_VariantReads,Description=\"Number of variant reads is less than 3\">"
            writer.puts "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than 0.1\">"
            writer.puts "##FILTER=<ID=single_strand,Description=\"All variant reads are in a single strand direction\">"
            writer.puts "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than 6 (default)\">"
            writer.puts "##FILTER=<ID=No_data,Description=\"No valid reads on this site\">"
            writer.puts "##FILTER=<ID=No_var,Description=\"No valid variant read on this site\">"
            writer.puts "##FILTER=<ID=ReqIncl,Description=\"Site does not pass filtering requirements, but is in the list of sites required to be included\">"
            writer.puts "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t#{@samples.keys.join("\t")}"
            STDOUT.puts "SAMPLENAMNE: #{@samples.keys.join("\t")}"
 
            #Printing column names
            @variants.each do |var, varObj|
              writer.puts "#{varObj.chr}\t#{varObj.coor}\t#{varObj.id}\t#{varObj.refBase}\t#{varObj.altBase}\t#{varObj.qual}\t#{varObj.filter}\t.\t#{varObj.format}"
            end
            writer.close
          end
          commandLineSort(out_file)
	end
	
	private

	def fileDoesNotCheckout(file)
          if File.exist?(file)
            return false
          else
            raise "Either #{file} does not exist."
            return true
          end
	end
	
	def commandLineSort(outFile)
          #Using combination if sed,grep and sort to sort the vcf file
          tfile = Tempfile.new('vcfPrinter')
          #substitute X,Y,M with 23,24,25 using sed
          `sed -i -e 's/^X/23/' -e 's/^Y/24/' -e 's/^M/25/' #{outFile}`
          #copy header
          `grep -e '^#' #{outFile} >#{tfile.path}`
          #sort rows
	  STDERR.puts 'Sorting sites file...'
          #`grep -v -e '^#' #{outFile} | sort -t  $'\t' -k1 -k2 -n >>#{tfile.path}`
          `grep -v -e '^#' #{outFile} | sort  -k1 -k2 -n >>#{tfile.path}`
          #revert back to X,Y,M from 23,24,25
          `sed -i -e 's/^23/X/' -e 's/^24/Y/' -e 's/^25/M/' #{tfile.path}`
          FileUtils.mv(tfile.path, outFile)
          STDERR.puts 'Sorting done.'
          return outFile
        end
end

if __FILE__== $0
	#puts Dir.glob(ARGV[0])
	makeSiteObj= MakeSiteFile.new(Dir.glob(ARGV[0]), ARGV[1], false ,'1:1-1000000')
end
