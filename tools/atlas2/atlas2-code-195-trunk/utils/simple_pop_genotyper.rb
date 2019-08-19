#!/usr/bin/ruby



class Simple_pop_genotyper

	def initialize( cutoffs=[0.05,0.9,0.9,0.2,2] )
		@heter_cutoff = cutoffs[0]
		@homo_var_cutoff = cutoffs[1]
		@homo_ref_cutoff = cutoffs[2]
		@ambig_cutoff = cutoffs[3]
		@depth_cutoff = cutoffs[4]
	end


	def genotype(ref_depth, alt_depths, total_depth, allele=1)
		major_allele = -1
		minor_allele = nil
		major_depth = ref_depth
		minor_depth = -1
		alt_depths.each_with_index do |alt, i|
			if(alt > major_depth)
				minor_allele = major_allele
				minor_depth = major_depth
				major_allele = i
				major_depth = alt
			elsif(alt > minor_depth)
				minor_allele = i
				minor_depth = alt
			end
		end
		major_allele += 1
		minor_allele += 1
		major_depth = major_depth.to_f
		minor_depth = minor_depth.to_f
		total_depth = total_depth.to_f
		
		if(major_allele == 0 || minor_allele == 0)
		   	allele = nil
			alt_depth = nil
			if major_allele != 0
				alt_depth = major_depth
				allele = major_allele
			end
			if minor_allele != 0
				alt_depth = minor_depth
				allele = minor_allele
			end


			if alt_depth == 0
				n = total_depth
			else
				alt_depth = alt_depth.to_f
				n = ref_depth + alt_depth
			#	return "#{allele}/." if total_depth < @depth_cutoff
			end
			return './.' if total_depth < @depth_cutoff
			if(n/total_depth < @ambig_cutoff)
				return './.'
			elsif(alt_depth/n >= @homo_var_cutoff && alt_depth >= @depth_cutoff)
				if(alt_depth < 2*@depth_cutoff)
					return "#{allele}/."
				else
					return "#{allele}/#{allele}"
				end
			elsif(alt_depth/n > @heter_cutoff && alt_depth >= @depth_cutoff)
				if(ref_depth/n > @heter_cutoff && ref_depth >= @depth_cutoff)
					return "#{allele}/0"
				else
					return "#{allele}/."
				end
			elsif(ref_depth >= @depth_cutoff && ref_depth/n >= @heter_cutoff)
				if(ref_depth/n >= @homo_ref_cutoff && ref_depth >= 2*@depth_cutoff)
					return '0/0'
				else
					return '0/.'
				end
			else
				return './.'
			end
		else
			if minor_depth == 0
				n = total_depth
			else
				n = major_depth + minor_depth
			end
			return './.' if total_depth < @depth_cutoff
			if(n/total_depth < @ambig_cutoff)
				return './.'
			elsif(major_depth/n >= 0.9 && major_depth >= @depth_cutoff)
				if(major_depth < 2*@depth_cutoff)
					return "#{major_allele}/."
				else
					return "#{major_allele}/#{major_allele}"
				end
			elsif(minor_depth/n > @heter_cutoff && minor_depth >= @depth_cutoff)
				return "#{major_allele}/#{minor_allele}"
			elsif(major_depth/n > @heter_cutoff && major_depth >= @depth_cutoff)
				return "#{major_allele}/."
			else
				return './.'
			end
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
	if(ARGV.length < 3)
		puts "USAGE: ruby simple_pop_genotyper.rb [vcf_file] coverage_file samples_list [-g (include genotype probabilities)]"
		exit(1)
	end

	vcf_file = ARGV[0]
	coverage_file = ARGV[1]
	samples_file = ARGV[2]
	include_quals = true if ARGV[3] == '-g'
	samples = Hash.new
	sample_names = Array.new
	File.open(samples_file, 'r').each_line do |line|
		sample = line.chomp
		samples[sample] = Hash.new
		sample_names.push sample
	end

	genotyper = Simple_pop_genotyper.new()
	labels=nil
	cov_reader = File.open(coverage_file, 'r')
	File.open(vcf_file, 'r').each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				puts "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
				puts '##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Reference Read Depth">'
				puts '##FORMAT=<ID=VR,Number=.,Type=Integer,Description="Variant Read Depth(s)">'
				labels = line.split("\t") 
				if(include_quals)
					puts "##FORMAT=<ID=PHR,Number=1,Type=Float,Description=\"Probability of a homozygous reference genotype\">"
					puts "##FORMAT=<ID=PH,Number=1,Type=Float,Description=\"Probability of a heterozygous genotype\">"
					puts "##FORMAT=<ID=PHV,Number=1,Type=Float,Description=\"Probability of a homozygous variant genotype\">"
				end
				puts "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	#{sample_names.join("\t")}"
				next
			end
			puts line
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)
		vcf_line.samples=samples
		vcf_line.sample_names = sample_names
		if(include_quals)
			vcf_line.format= [:"GT",:"PHR",:"PH",:"PHV",:"VR",:"RR",:"DP"]
		else
			vcf_line.format= [:"GT",:"VR",:"RR",:"DP"]
		end
		cov_line = cov_reader.gets
		cov_line.chomp!
		covs = cov_line.split(',')
		raise "incorrect coverage line position #{cov_line}\n#{vcf_line.to_s}" if covs[0] != vcf_line.chr.to_s || covs[1].to_i != vcf_line.pos
		i = 2
		alt_ids = Array.new
		vcf_line.alt.each do |alt|
			ref_s = vcf_line.ref.to_s
			alt_s = alt.to_s
			if(ref_s.length > alt_s.length) # deletion
				alt_ids.push("[#{(alt_s.length - ref_s.length).to_s}]")
			else # insersion
				l = alt_s.length - ref_s.length
				alt_ids.push("[#{alt_s[1..l]}]")
			end
		end
		vcf_line.samples.each_value do |sample_data|
			alt_depths = Array.new
			raise "parser got lost! possibly missing cdb file #{covs[i]}" if ! (covs[i] =~ /S\d/)
			i+=1
			total_depth = covs[i].to_i
			i+=1
			ref_depth = covs[i].to_i
			i+=1
			while( i< covs.length && !(covs[i] =~ /S\d/) )
				allele_i = alt_ids.index(covs[i])
				raise "could not find alt #{covs[i]} in the coverage file!\n #{vcf_line.print}" if allele_i.nil?
				i+=1
				alt_depths.push(covs[i].to_i)
				i+= 1
			end
		
			# replace '.' with 0
			sample_data[:"VR"]=alt_depths.join(',')
			#sample_data[:"VR"]=alt_depths.max
			sample_data[:"RR"]=ref_depth
			sample_data[:"DP"]=total_depth
			geno = genotyper.genotype(ref_depth, alt_depths, total_depth)
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
