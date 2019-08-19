$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 2)
	puts "USAGE: ruby add_source_qaul.rb vcf_file source_vcf"
	exit 0
end

source_vcf = ARGV[1]
labels=nil
prev_line = nil	
header_line = '##INFO=<ID=SFQ,Number=.,Type=String,Description="Source file variant quality score">'
has_header = false
source_header = `grep "#CHROM	" #{source_vcf}`
File.open(ARGV[0], 'r').each_line do |line|
	line.chomp!
	if(line[0,1] == '#')
		has_header = true if line == header_line
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
			puts header_line unless has_header
		end
		puts line
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	source_str = `grep "#{vcf_line.chr}	#{vcf_line.pos}	" #{source_vcf}`
	source_str.chomp!
	qual = 0
	if( (!source_str.nil?) && source_str.length > 0 ) # this line is in the source VCF
		cols = source_str.split("\t")
		qual = cols[5].to_i
	end
	if( qual == 0 && !(line.include?('illum')) )
		qual = 'NA'
	end
	if(vcf_line.info[:"SFQ"].nil?)
		vcf_line.info[:"SFQ"] = qual
	elsif(vcf_line.info[:"SFQ"].class != Array)
		tmp = vcf_line.info[:"SFQ"]
		vcf_line.info[:"SFQ"] = [tmp,qual]
	elsif(vcf_line.info[:"SFQ"].class == Array)
		vcf_line.info[:"SFQ"].push(qual)
	else
		raise "unhandled case, unknown SF entry: #{vcf_line.info[:"SF"]}"
	end
	puts vcf_line.print
end
