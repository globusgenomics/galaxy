#!/usr/bin/ruby


class Simple_genotyper

	def initialize( cutoffs=[0.1,0.9,0.9,0.5,4] )
		@heter_cutoff = cutoffs[0]
		@homo_var_cutoff = cutoffs[1]
		@homo_ref_cutoff = cutoffs[2]
		@ambig_cutoff = cutoffs[3]
		@depth_cutoff = cutoffs[4]
	end


	def genotype(ref_depth, alt_depth, total_depth, allele=1)
		ref_depth = ref_depth.to_f
		total_depth = total_depth.to_f
		
		if alt_depth == 0
			n = total_depth
			return './.' if total_depth < @depth_cutoff
		else
			alt_depth = alt_depth.to_f
			n = ref_depth + alt_depth
			return "#{allele}/." if total_depth < @depth_cutoff
		end
		if(n/total_depth < @ambig_cutoff)
			return './.'
		elsif(alt_depth/n >= @homo_var_cutoff)
			return "#{allele}/#{allele}"
		elsif(alt_depth/n > @heter_cutoff)
			return './.' if allele == '.'
			return "#{allele}/0"
		elsif(ref_depth/n >= @homo_ref_cutoff)
			return '0/0'
		else
			return './.'
		end
	end
end

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
	if(ARGV.length < 1)
		puts "USAGE: ruby simple_genotyper.rb [vcf_file] [-g (include genotype probabilities)]"
		exit(0)
	end

	vcf_file = ARGV[0]
	include_quals = true if ARGV[1] == '-g'

	genotyper = Simple_genotyper.new()
	labels=nil
	File.open(vcf_file, 'r').each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
				if(include_quals)
					puts "##FORMAT=<ID=PHR,Number=1,Type=Float,Description=\"Probability of a homozygous reference genotype\">"
					puts "##FORMAT=<ID=PH,Number=1,Type=Float,Description=\"Probability of a heterozygous genotype\">"
					puts "##FORMAT=<ID=PHV,Number=1,Type=Float,Description=\"Probability of a homozygous variant genotype\">"
				end
			end
			puts line
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)
		if(include_quals)
			vcf_line.format= [:"GT",:"PHR",:"PH",:"PHV",:"VR",:"RR",:"DP"]
		end
		vcf_line.samples.each_value do |sample_data|
			alt_depth = get_depth(sample_data[:"VR"])
			ref_depth = get_depth(sample_data[:"RR"])
			total_depth = get_depth(sample_data[:"DP"])
			# replace '.' with 0
			sample_data[:"VR"]=alt_depth
			sample_data[:"RR"]=ref_depth
			sample_data[:"DP"]=total_depth
			old_geno = sample_data[:"GT"]
			if(old_geno =~ /[123456789]/)
				allele = $&
			else
				allele = '.'
			end
			geno = genotyper.genotype(ref_depth, alt_depth, total_depth, allele)
			sample_data[:"GT"]=geno 
			if(include_quals)
				if(total_depth > 3)
					if(geno == '0/0')
						sample_data[:"PHR"]=0.95
						sample_data[:"PH"]=0.025
						sample_data[:"PHV"]=0.025
					elsif(geno =~ /[123456789]\/[123456789]/)
						sample_data[:"PHR"]=0.025
						sample_data[:"PH"]=0.025
						sample_data[:"PHV"]=0.95
					elsif(geno =~ /[123456789]\/0/)
						sample_data[:"PHR"]=0.025
						sample_data[:"PH"]=0.95
						sample_data[:"PHV"]=0.025
					else # ambiguous genotype/no data
						sample_data[:"PHR"]=0.34
						sample_data[:"PH"]=0.33
						sample_data[:"PHV"]=0.33
					end
				else
					if(geno =~ /[123456789]\/\./)
						sample_data[:"PHR"]=0.025
						sample_data[:"PH"]=0.4875
						sample_data[:"PHV"]=0.4875
					else # genotype must be ./.
						sample_data[:"PHR"]=0.34
						sample_data[:"PH"]=0.33
						sample_data[:"PHV"]=0.33
					end

				end
			end
		end
		puts vcf_line.print()
	end
end
