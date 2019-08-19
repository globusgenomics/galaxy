#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Variant_af_by_sample

	def initialize(file)
		@filename = file
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
			process_line(vcf_line)
		end
	end

	def process_line(line)
		line.sample_names.each_with_index do |sample_name, sample_index|
			return if line.filter.to_s != "PASS"
			puts line.info[:"AAF"] if line.samples[sample_name][:"GT"].to_s =~ /[1-9]/
		end
	end

end

if(ARGV.length < 1)
	puts "USAGE: ruby variant_af_by_sample.rb vcf_file"
	exit 1
end

extractor = Variant_af_by_sample.new(ARGV[0])
extractor.run
