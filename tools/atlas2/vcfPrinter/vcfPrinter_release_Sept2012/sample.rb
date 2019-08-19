require 'pileup.rb'
require 'vcf.rb'

class Sample
	
	attr_reader :name, :pileup, :vcf

	def initialize(name, vcfFile)
		@name=name
		#@pileup = Pileup.new(bamFile, bedFile, reference)
		@vcf = Vcf.new(vcfFile)
		@pileupExist=false
	end

	def generatePileup(bamFile,bedFile,reference)
		puts "#{bamFile}\n#{bedFile}\n#{reference}"
		@pileup = Pileup.new(false,bamFile, bedFile, reference)
		@pileupExist=true
	end

	def createPileup(pileupFile)
		@pileup = Pileup.new(true,pileupFile)
	end

	def hasPileup
		return @pileupExist
	end

end

