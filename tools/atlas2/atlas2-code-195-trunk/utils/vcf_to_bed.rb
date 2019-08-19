$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 1)
	puts "USAGE: ruby vcf_to_bed.rb <vcf_file>"
	exit 0
end

labels=nil
prev_line = nil	
File.open(ARGV[0], 'r').each_line do |line|
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
		end
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	puts "#{vcf_line.chr} #{vcf_line.pos} #{vcf_line.pos+1}"
end
