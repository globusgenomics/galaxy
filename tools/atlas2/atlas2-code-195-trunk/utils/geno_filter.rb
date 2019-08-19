#!/usr/bin/ruby


if __FILE__ == $0

	def get_depth(entry)
		if( entry == '.')
			return 0
		else
			return entry.to_i
		end
	end

	$:.unshift File.join(File.dirname(__FILE__),'.')
	require 'vcf_line'
	if(ARGV.length < 2)
		puts "USAGE: cat [vcf-file] | ruby geno_filter.rb [geno qual cut-off] [depth cut-off]"
		exit(0)
	end

	#vcf_file = ARGV[0]
	gq_cutoff = ARGV[0].to_i
	depth_cutoff = ARGV[1].to_i

	labels=nil
	STDIN.each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
			end
			puts line
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)
		#vcf_line.format.push(:"FT") unless vcf_line.format.include?(:"FT")
		vcf_line.info[:"NS"] = 0  #initialize NS counts
		vcf_line.samples.each_value do |sample_data|
			gq = sample_data[:"GQ"]
			if(gq.nil? || gq == '.')
		#		sample_data[:"FT"]='.'
				next
			elsif(gq.to_i < gq_cutoff)
				sample_data[:"GT"] = './.'
		#		sample_data[:"FT"]='genoQual'
			elsif(sample_data[:"DP"].to_i < depth_cutoff)
				sample_data[:"GT"] = './.'
		#		sample_data[:"FT"]='minDepth'
			else
				vcf_line.info[:"NS"] += 1 #recount NS
		#		sample_data[:"FT"]='PASS'
			end
		end
		puts vcf_line.print()
#		puts vcf_line.info[:"NS"]
	end
end
