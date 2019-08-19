#!/usr/bin/ruby


if __FILE__ == $0


	$:.unshift File.join(File.dirname(__FILE__),'.')
	require 'vcf_line'
	if(ARGV.length < 1)
		puts "USAGE: ruby indel_metrics.rb vcf-file"
		exit(0)
	end

	vcf_file = ARGV[0]

	labels=nil
	var_sample_site_count = 0
	indel_count = 0
	frameshift_count = 0
	total_alleles_count = 0
	frameshift_allele_count = 0
	homo_fs_count = 0
	heter_fs_count = 0 
	homo_inframe_count = 0
	heter_inframe_count = 0
	one_bp_count = 0.0
	two_bp_count = 0.0
	three_bp_count = 0.0
	ins_count = 0
	del_count = 0
	prev_locus = nil
	`rm -f #{ARGV[0]}.allele_lengths`
	`rm  -f #{ARGV[0]}.lengths`
	File.open(vcf_file, 'r').each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
			end
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)
		allele_fs = Array.new
		indel_count += 1
		vcf_line.alt.each do |allele|
			if ( ( allele.to_s.length - vcf_line.ref.to_s.length ) % 3 == 0 )
				allele_fs.push(false)
			else
				allele_fs.push(true)
				frameshift_allele_count += 1
			end
			total_alleles_count += 1
			length = allele.to_s.length - vcf_line.ref.to_s.length
			`echo #{length} >> #{ARGV[0]}.allele_lengths`
			if(length > 0)
				ins_count += 1
			elsif(length < 0)
				del_count += 1
			else
				raise "not an indel! \n#{line.to_s}"
			end
		end
		vcf_line.samples.each do |name, sample_data|
			geno = sample_data[:"GT"]
			if(geno.nil? || geno == './.' || geno == '0/0' || geno == '0/.' || geno == './0' || geno == '.' || geno == '0' || geno == '.|.' || geno == '0|0' || geno == "0|." || geno == ".|0")
				next
			end
			var_sample_site_count += 1

			# get the right allele for the sample
			if(geno =~ /[123456789]/)
				allele = $&
			else 
				raise "could not find allele number in geno: #{geno}"
			end
			length = vcf_line.alt[allele.to_i-1].to_s.length - vcf_line.ref.to_s.length
			#`echo #{length} >> #{ARGV[0]}.lengths`
			if(geno == "#{allele}/#{allele}" || geno == "#{allele}|#{allele}") 
				is_homo = true
			else
				is_homo = false
			end
			if( length.abs == 1 )
				if(is_homo)
					one_bp_count += 1
				else
					one_bp_count += 1
				end
			end
			if( length.abs == 2 )
				if(is_homo)
					two_bp_count += 1
				else
					two_bp_count += 1
				end
			end
			if( length.abs == 3 )
				if(is_homo)
					three_bp_count += 1
				else
					three_bp_count += 1
				end
			end
			if( allele_fs[(allele.to_i)-1] ) # allele is a frameshift
				frameshift_count += 1
				if(is_homo)
					homo_fs_count += 1
				else
					heter_fs_count += 1
				end
			else
				if(is_homo)
					homo_inframe_count += 1
				else
					heter_inframe_count += 1
				end
			end
		end
		prev_locus = vcf_line.locus
	end
	sample_count = labels.length - 9
	homo_fs_rate = homo_fs_count.to_f/(homo_fs_count + heter_fs_count).to_f
	homo_inframe_rate = homo_inframe_count.to_f/(homo_inframe_count + heter_inframe_count).to_f
	puts "# Indels	Indels-per-sample	Total SS	FS_rate_by_allele	FS_rate_mean_per_indiv	Homo FS rate	Homo inframe rate	Homo inframe/FS ratio	1bp/2bp Allele Ratio	3bp/2bp Allele Ratio	1bp/3bp_allele_ratio	num_dels	num_ins	del_ins_ratio"
	puts "#{indel_count}	#{var_sample_site_count.to_f/sample_count.to_f}	#{var_sample_site_count}	#{frameshift_allele_count.to_f/total_alleles_count.to_f}	#{frameshift_count.to_f/var_sample_site_count.to_f}	#{homo_fs_rate}	#{homo_inframe_rate}	#{homo_inframe_rate/homo_fs_rate}	#{one_bp_count/two_bp_count}	#{three_bp_count/two_bp_count}	#{one_bp_count/three_bp_count}	#{del_count}	#{ins_count}	#{del_count.to_f/ins_count.to_f}"
end
