#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Sample_site_printer

	def initialize(file, pos, samples)
		@filename = file
		@chr,@pos= pos.split(":")
		@samples=samples.split(",")
	end

	def run
		labels=nil
		File.open(@filename, 'r').each_line do |line|
			line.chomp!
			if(line[0,1] == '#')
				if( line[0,4] == '#CHR' )
					labels = line.split("\t") 
				end
				next
			end
			vcf_line = Vcf_line.read_line(line, false, labels)
			if(vcf_line.chr.to_s == @chr && vcf_line.pos == @pos.to_i)
				print "#{vcf_line.chr}	#{vcf_line.pos}	#{vcf_line.ref}	#{vcf_line.alt.join(',')}	AC=#{vcf_line.info[:"AC"]}"
				@samples.each do |sample|
					print "	#{sample}:#{vcf_line.print_sample(sample)}"
				end
				print "\n"
				break
			end
		end
	end

end

if(ARGV.length<3)
	puts "USAGE: print_sample_sites.rb vcf chr:pos sample(s)"
	exit 1
end
printer = Sample_site_printer.new(ARGV[0], ARGV[1], ARGV[2])
printer.run



