$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 1)
	puts "USAGE: ruby merge_duplicate_alts vcf_file"
	exit 0
end

labels=nil
File.open(ARGV[0], 'r').each_line do |line|
	line.chomp!
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
		end
		puts line
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	first_ref_base = vcf_line.ref.to_s[0]
	alts = vcf_line.alt
	alts.each_with_index do |alt, i|
		alt_str = alt.to_s
		if(alt_str[0] != first_ref_base)
			alt_str[0]=first_ref_base
			alts[i]=alt_str.to_sym
		end
	end
	puts vcf_line.print
end
