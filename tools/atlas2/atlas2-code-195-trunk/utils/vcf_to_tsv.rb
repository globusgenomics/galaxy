$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'


def get_word(ref, alt)
	new_ref, new_alt = Vcf_line.simplify(ref.to_s, alt.to_s)
	if(new_alt.length > new_ref.length) # insertion
		new_alt[0]=""
		return new_alt
	else
		new_ref[0]=""
		return new_ref
	end
end

if(ARGV.length < 1)
	puts "USAGE: ruby vcf_to_tsv.rb <vcf_file> <comma-delimited INFO IDs>"
	exit 0
end

if(ARGV[1].nil?)
	info_fields = []
else
	info_fields = ARGV[1].split(',')
end
labels=nil
prev_line = nil	
print "chr	coor	ref	alt	seq	length	qual	filter	"
puts info_fields.join("\t")
File.open(ARGV[0], 'r').each_line do |line|
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
		end
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	if(vcf_line.alt.length > 1)
		0.upto(vcf_line.alt.length-1) do |allele_i|
			print "#{vcf_line.chr}	#{vcf_line.pos}	#{vcf_line.ref}	#{vcf_line.alt[allele_i]}	#{get_word(vcf_line.ref, vcf_line.alt[allele_i])}	#{vcf_line.alt[allele_i].to_s.length-vcf_line.ref.to_s.length}	#{vcf_line.qual}	#{vcf_line.filter}"
			info_fields.each do |id|
				print "	#{vcf_line.info[id.to_sym][allele_i]}"
			end
			print "\n"
		end
	else
		print "#{vcf_line.chr}	#{vcf_line.pos}	#{vcf_line.ref}	#{vcf_line.alt[0]}	#{get_word(vcf_line.ref, vcf_line.alt[0])}	#{vcf_line.alt[0].to_s.length-vcf_line.ref.to_s.length}	#{vcf_line.qual}	#{vcf_line.filter}"
		info_fields.each do |id|
			print "	#{vcf_line.info[id.to_sym]}"
		end
		print "\n"
	end
end
