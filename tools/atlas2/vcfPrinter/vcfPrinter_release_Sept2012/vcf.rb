require 'vcf_variant.rb'

class Vcf

	attr_accessor :vcf

	def initialize(vcfFile)
		@vcf = Hash.new
		File.open(vcfFile, 'r').each_line do |line|
		if line =~ /^#/
 	                next
                else
			vcfObj = VcfLine.new(line)
			@vcf[vcfObj.uniqPos]=vcfObj
		end
		end
	end
end


