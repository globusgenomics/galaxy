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
		puts "USAGE: remove_non-snps.rb <vcf file>\n note: does not work on sites only files"
		exit(0)
	end

	vcf_file = ARGV[0]

	labels=nil
	File.open(vcf_file, 'r').each_line do |line|
		skip = false
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
				puts "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Allele count in genotypes\">"
				puts "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">"
				puts "##INFO=<ID=AAF,Number=.,Type=Float,Description=\"Alternative allele frequency\">"
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
		vcf_line.sample_names.each do |sample_name|
			geno = vcf_line.samples[sample_name][:"GT"]
			if(geno =~ /([\d\.]).([\d\.])/)
				a1 = $1
				a2 = $2
				if(a1 =~ /\d/)
					if(vcf_line.alt[a1.to_i-1].to_s.length == vcf_line.ref.to_s.length)
						an += 1
						if(a1 != '0')
							if ac[a1.to_i-1] == nil
								STDERR.puts "WARNING: could not find allele #{a1} in sample #{sample_name} on vcf line: #{vcf_line.chr}:#{vcf_line.pos}   skipping line..."
								skip = true
								break
							end
							ac[a1.to_i-1] += 1
						end
					else
						vcf_line.samples[sample_name][:"GT"].sub!(a1, '.')
					end
				end
				if(a2 =~ /\d/)
					if(vcf_line.alt[a2.to_i-1].to_s.length == vcf_line.ref.to_s.length)
						an += 1
						if(a2 != '0')
							if ac[a2.to_i-1] == nil
								STDERR.puts "WARNING: could not find allele #{a2} in sample #{sample_name} on vcf line: #{vcf_line.chr}:#{vcf_line.pos}   skipping line..."
								skip = true
								break
							end
							ac[a2.to_i-1] += 1
						end
					else
						vcf_line.samples[sample_name][:"GT"].sub!(a2, '.')
					end
				end
			end
		end
		next if skip
		
		if(!vcf_line.sample_names.nil? && vcf_line.sample_names.length > 0)
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
		vcf_line.info[:"AAF"] = mafs
		if(ac.length == 0)
			#STDERR.puts "Line #{vcf_line.chr}:#{vcf_line.pos} has no variant calls, it will be removed."
		else
			puts vcf_line.print()
		end
	end
end
