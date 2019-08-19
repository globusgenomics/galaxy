#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Variant_sample_printer

	def initialize(file, max_samples)
		@filename = file
		@max_samples = max_samples.to_i
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
			i=0
			print "#{vcf_line.chr}	#{vcf_line.pos}	#{vcf_line.ref}	#{vcf_line.alt.join(',')}	AC=#{vcf_line.info[:"AC"]}	"
			vcf_line.sample_names.each do |name|
				break if(i>= @max_samples)
				if(vcf_line.samples[name][:"GT"].to_s =~ /[1-9]/)
					print ',' unless i==0
					print name
					i+=1
				end
			end	
			print "\n"
		end
	end

end

if(ARGV.length<2)
	puts "USAGE: print_variant_samples.rb vcf max_samples"
	exit 1
end
printer = Variant_sample_printer.new(ARGV[0], ARGV[1])
printer.run



