#!/usr/bin/ruby


if __FILE__ == $0
	$:.unshift File.join(File.dirname(__FILE__),'.')
	require 'vcf_line'
	require 'ref_seq'

	if(ARGV.length < 2)
		puts "USAGE: ruby homopolymer_filter.rb vcf_file reference_genome"
		exit(0)
	end

	vcf_file = ARGV[0]
	ref_genome = ARGV[1]
	reference = Ref_seq.new(ref_genome)
	labels=nil
	File.open(vcf_file, 'r').each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
			end
			puts line
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)
		#next if vcf_line.alt.length > 1 # skip multi-allelic sites
		vcf_line.alt.each do |allele|
			ref = vcf_line.ref.to_s
			alt= allele.to_s
			ref,alt = Vcf_line.simplify(vcf_line.ref.to_s,allele.to_s) if ref.length != 1 && alt.length != 1
			substring = ""
			if(alt.length < ref.length) # deletion
				substring = ref
				substring += reference.get_base(vcf_line.chr, vcf_line.pos+ref.length)
				substring.squeeze!
				if substring.length > 2
					puts line
					break
				end
			else  # insertion
				substring = alt
				substring += reference.get_base(vcf_line.chr, vcf_line.pos+1)
				substring.squeeze!
				if substring.length > 2
					puts line
					break
				end

			end
		end
	end
end
