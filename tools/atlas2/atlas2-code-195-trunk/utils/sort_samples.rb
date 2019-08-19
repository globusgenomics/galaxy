$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 1)
	puts "USAGE: ruby sort_samples.rb vcf_line"
	exit 0
end

window = ARGV[1].to_i
labels=nil
prev_line = nil	
sample_names = nil
File.open(ARGV[0], 'r').each_line do |line|
	line.chomp!
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
			basic_labels = labels[0..8]
			sample_names = labels[9...labels.length]
			sample_names.sort!
			puts "#{basic_labels.join("\t")}	#{sample_names.join("\t")}"
		else
			puts line
		end
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	vcf_line.sample_names = sample_names
	puts vcf_line.print
end
