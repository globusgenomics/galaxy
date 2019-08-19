require 'pileup_variant.rb'


class Pileup

	attr_accessor :pileup

	def initialize(hasPileup=false, *args)
		@pileup = Hash.new()
		if hasPileup
			output = File.open(args[0],'r')
			output.each_line do |out|
			pileupObj = PileupLine.new(out)
			@pileup[pileupObj.pos]=pileupObj
			end
		else
			#puts "samtools pileup -l #{args[1]} -f #{args[2]} #{args[0]}"
			puts "samtools mpileup -l #{args[1]} -f #{args[2]} #{args[0]}"
			#`samtools pileup -l #{args[1]} -f #{args[2]} #{args[0]}`.each_line do |out|  #modified by liubo
			`samtools mpileup -l #{args[1]} -f #{args[2]} #{args[0]}`.each_line do |out|
			pileupObj = PileupLine.new(out)
			@pileup[pileupObj.pos]=pileupObj
			end
		end
    	end
end

