#!/usr/bin/ruby

def get_repeat_substr(seq)
	if(seq =~ /(.+?)\1+/) # find shortest repeat unit
		if($&.length == seq.length)
			return $1
		end
	end
	return seq
end

if __FILE__ == $0
	$:.unshift File.join(File.dirname(__FILE__),'.')
	require 'vcf_line'

	if(ARGV.length < 1)
		puts "USAGE: ruby mark_frameshift.rb vcf_file"
		exit(0)
	end

	vcf_file = ARGV[0]
	labels=nil
	File.open(vcf_file, 'r').each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
				puts "##INFO=<ID=FS,Number=A,Type=Character,Description=\"Does the allele cause a frameshift? (T/F)\">"
				puts "##INFO=<ID=FS,Number=A,Type=Integer,Description=\"Variant length\">"
			end
			puts line
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)
		vcf_line.info[:"FS"]=Array.new
		vcf_line.info[:"length"]=Array.new
		vcf_line.alt.each do |allele|
			ref = vcf_line.ref.to_s
			alt= allele.to_s
			ref,alt = Vcf_line.simplify(vcf_line.ref.to_s,allele.to_s) if ref.length != 1 && alt.length != 1
			length = alt.length - ref.length 
			vcf_line.info[:"length"].push(length)
			if( length % 3 == 0)
				vcf_line.info[:"FS"].push("F")
			else
				vcf_line.info[:"FS"].push("T")
			end
		end
		puts vcf_line.print
	end
end
