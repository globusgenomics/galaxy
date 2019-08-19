
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
require 'bed_file.rb'

if(ARGV.length < 2)
	puts "USAGE: ruby bed_filter.rb vcf_file bed_file"
	exit 0
end

labels=nil
bed = Bed_file.new(ARGV[1])
	
File.open(ARGV[0], 'r').each_line do |line|
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
		end
		puts line
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	chr = vcf_line.chr
	chr.slice!('chr')
	if(bed.pos_included?(chr.to_sym, vcf_line.pos))
		puts line
	end
end
		

