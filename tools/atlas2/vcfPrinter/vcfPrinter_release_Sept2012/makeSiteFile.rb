require 'vcf_variant.rb'
require 'tempfile'

$:.unshift File.join(File.dirname(__FILE__),'.')


class MakeSiteFile
	attr_reader :samples, :variants

	def initialize(vcf_files, out_file, onlyPASS=false)
		@samples = Hash.new
                @variants = Hash.new
                vcf_files.each do |vcf_file|
			next if fileDoesNotCheckout(vcf_file)
			#puts vcf_file
			@samples[File.basename(vcf_file, '.vcf')]=vcf_file
			#puts vcf_file
		end
		puts "Number of Samples: #{@samples.length}"
                @samples.each do |sampleName, fileName|
			lines = File.open(fileName, 'r')
                        lines.each do |line|
                        next if line =~ /^#/
                        var = VcfLine.new(line)
			if onlyPASS
				next if var.filter!='PASS'
				@variants[var.uniqPos].updateQual(var.qual) if !@variants[var.uniqPos].nil?
				@variants[var.uniqPos]=var
			else
				@variants[var.uniqPos].updateQual(var.qual) if !@variants[var.uniqPos].nil?
				@variants[var.uniqPos]=var
			end
			end
			lines.close
		end
		puts "Number of Variants: #{@variants.length}"
                File.open(out_file,"w") do |writer|
                #Printing VCF header information
                        writer.puts "##fileformat=VCFv4.0"
                        time = Time.new
                        writer.puts "##fileDate=#{time.ctime}"
                        writer.puts "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
                        writer.puts "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
                        writer.puts "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
                        writer.puts "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Variant Read Depth\">"
                        writer.puts "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">"
                        writer.puts "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">"
                        writer.puts "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
                        writer.puts '##FORMAT=<ID=FT,Number=1,Type=String,Description="Sample Genotype Filter">'
                        writer.puts "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t#{@samples.keys.join("\t")}"
			#Printing column names
			@variants.each do |var, varObj|
                                writer.puts "#{varObj.chr}\t#{varObj.coor}\t#{varObj.id}\t#{varObj.refBase}\t#{varObj.altBase}\t#{varObj.qual}\t#{varObj.filter}\t.\tGT:VR:RR:DP:GQ:FT"
                        end
		writer.close
		end
		commandLineSort(out_file)
	end
	
	private

	def fileDoesNotCheckout(file)
		if File.exist?(file) and File.fnmatch("*.vcf",file)
			return false
                else
                	raise "Either #{file} does not exist or the file does not have proper .vcf extension"
			return true
			#return false
                end
	end
	
	def commandLineSort(outFile)
                tfile = Tempfile.new('tmp')
                #substitute X,Y,M with 23,24,25 using sed
                `sed -i -e 's/^X/23/' -e 's/^Y/24/' -e 's/^M/25/' #{outFile}`
                #copy header
                `grep -e '^#' #{outFile} >#{tfile.path}`
                #sort rows
		puts 'Sorting sites file...'
#                `sleep 3600`
                #`grep -v -e '^#' #{outFile} | sort -t  $'\t' -k1 -k2 -n >>#{tfile.path}`
                `grep -v -e '^#' #{outFile} | sort  -k1 -k2 -n >>#{tfile.path}`
                #revert back to X,Y,M from 23,24,25
                `sed -i -e 's/^23/X/' -e 's/^24/Y/' -e 's/^25/M/' #{tfile.path}`
                FileUtils.mv(tfile.path, outFile)
                return outFile
        end
end

if __FILE__== $0
	#puts Dir.glob(ARGV[0])
	makeSiteObj= MakeSiteFile.new(Dir.glob(ARGV[0]), ARGV[1])
end
