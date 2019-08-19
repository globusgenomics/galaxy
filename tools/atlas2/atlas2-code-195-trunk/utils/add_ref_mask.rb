#!/usr/bin/ruby
$:.unshift File.join(File.dirname(__FILE__))
require 'vcf_line.rb'
require 'ref_seq.rb'


class Vcf_editor

	def initialize(file, mask_file, ascii)
		@filename = file
		@mask = Ref_seq.new(mask_file)
		@ascii = ascii
	end

	def run
		labels=nil
		File.open(@filename, 'r').each_line do |line|
			line.chomp!
			if(line[0,1] == '#')
				if( line[0,4] == '#CHR' )
					labels = line.split("\t") 
				end
				puts line
				next
			end
			vcf_line = Vcf_line.read_line(line, false, labels)
			edit_line(vcf_line)
			puts vcf_line.print()
		end
	end

	def edit_line(line)
		if(@ascii)
			qual = @mask.get_base(line.chr().to_sym,line.pos()).ord-33
			line.info()[:"Mask"]=qual
		else
			line.info()[:"Mask"]=@mask.get_base(line.chr().to_sym,line.pos())
		end
	end
end

if(ARGV.length < 1)
	puts "USAGE: ruby add_ref_mask.rb vcf_file ref_mask.fasta [-a (for acii QUAL score)]"
	exit 1
end

editor = Vcf_editor.new(ARGV[0], ARGV[1], ARGV[2]=="-a")
editor.run


