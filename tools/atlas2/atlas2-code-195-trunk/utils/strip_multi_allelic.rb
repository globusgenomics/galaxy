#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Multi_allele_stripper

	def initialize(file)
		@filename = file
		@id=0
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

	def edit_line(vcf_line)
		vcf_line.sample_names.each do |sample_name|
			vcf_line.samples[sample_name][:"GT"]=vcf_line.samples[sample_name][:"GT"].to_s.gsub(/[2-9]/, "1")
		end
		vcf_line.alt=[vcf_line.alt[0]]
	end
	
end

stripper = Multi_allele_stripper.new(ARGV[0])
stripper.run
