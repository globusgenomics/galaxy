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
	if(ARGV.length < 1)
		puts "USAGE: add_allele_freq.rb [vcf_file] -clean\n	the -clean option will remove alleles with an allele count (AC) of zero"
		exit(0)
	end

	vcf_file = ARGV[0]
	clean_alleles = (ARGV[1] == '-clean')

	labels=nil
	File.open(vcf_file, 'r').each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
				puts "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes\">"
				puts "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">"
				puts "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele frequency\">"
			end
			puts line
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)

		ac =  Array.new
		an = 0
		vcf_line.alt.each do |allele|
			ac.push 0
		end
		if(vcf_line.sample_names.length == 0)
			STDERR.puts "File #{ARGV[0]} has sites only, no samples.  No output produced."
			exit(1)
		end
		vcf_line.sample_names.each do |sample_name|
			unless(vcf_line.samples[sample_name].length == 0)
				# deal with sample filters
				if(vcf_line.samples[sample_name][:"FT"] != nil && vcf_line.samples[sample_name][:"FT"].to_s != "PASS" && vcf_line.samples[sample_name][:"FT"].to_s != ".")
					if(vcf_line.samples[sample_name][:"GT"] =~ /[1-9]/)
						vcf_line.samples[sample_name][:"GT"] = "./."
					elsif(vcf_line.samples[sample_name][:"FT"].to_s != "No_data") # don't list filters for non-variant genotypes (except No_data)
						vcf_line.samples[sample_name][:"FT"] = '.'
					end
				end

				geno = vcf_line.samples[sample_name][:"GT"]
				if(geno =~ /([\d\.]+)[\/\\|]([\d\.]+)/)
					a1 = $1
					a2 = $2
					if(a1 =~ /\d+/)
						an += 1
						if(a1 != '0')
							raise "could not find allele #{a2} in sample #{sample_name} on vcf line:\n #{line}" if ac[a1.to_i-1] == nil
							ac[a1.to_i-1] += 1
						end
					end
					if(a2 =~ /\d+/)
						an += 1
						if(a2 != '0')
							raise "could not find allele #{a2} in sample #{sample_name} on vcf line:\n #{line}" if ac[a2.to_i-1] == nil
							ac[a2.to_i-1] += 1
						end
					end
				end
			end
		end
		
		if(clean_alleles)
			new_ac = ac.dup
			removed_alleles = 0
			ac.each_with_index do |allele_count, i|
				if(allele_count == 0)
					vcf_line.alt.delete_at(i-removed_alleles)
					new_ac.delete_at(i-removed_alleles)
					removed_alleles += 1
				end	
				if(removed_alleles > 0)
					vcf_line.sample_names.each do |sample_name|
						allele_index = i + 1
						vcf_line.samples[sample_name][:"GT"].gsub!(allele_index.to_s, (allele_index - removed_alleles).to_s) unless vcf_line.samples[sample_name][:"GT"].nil?
					end
				end
			end
			ac = new_ac
		end
			
		mafs = Array.new
		vcf_line.info[:"AC"]=ac.join(",")
		vcf_line.info[:"AN"]=an
		ac.each do |one_ac|
			#if(one_ac.to_f/an > 0.5)
			#	mafs.push (1.0-one_ac.to_f/an)
			#else
			#	mafs.push (one_ac.to_f/an)
			#end
			mafs.push (one_ac.to_f/an)
		end	
		vcf_line.info[:"AF"] = mafs
		if(clean_alleles && ac.length == 0)
			STDERR.puts "Line #{vcf_line.chr}:#{vcf_line.pos} has no variant calls, it will be removed."
		else
			raise "Line #{vcf_line.chr}:#{vcf_line.pos} has no variant calls!" if ac.length == 0
			puts vcf_line.print()
		end
	end
end
