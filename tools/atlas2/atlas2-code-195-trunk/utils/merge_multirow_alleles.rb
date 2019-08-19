#!/usr/bin/ruby
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 1)
	puts "USAGE: ruby merge_multirow_alleles.rb [vcf_file]"
	exit 0
end

labels=nil
prev_line = nil	
merge_count = 0
File.open(ARGV[0], 'r').each_line do |line|
	line.chomp!
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
		end
		puts line
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	if(prev_line)
		chr = vcf_line.chr
		pos = vcf_line.pos
		if(chr == prev_line.chr && prev_line.pos == pos)
			merge_count += 1
			if(vcf_line.ref.to_s.length > prev_line.ref.to_s.length)
				allele_addition = vcf_line.ref.to_s[(prev_line.ref.to_s.length)...(vcf_line.ref.to_s.length)]
				prev_line.alt.each_with_index do |allele, i|
					prev_line.alt[i] = allele.to_s
					prev_line.alt[i] += allele_addition
				end
				prev_line.ref = vcf_line.ref
			else # prev line has larger ref
				allele_addition =  prev_line.ref.to_s[(vcf_line.ref.to_s.length)...(prev_line.ref.to_s.length)]
				vcf_line.alt.each_with_index do |allele, i|
					vcf_line.alt[i] = allele.to_s
					vcf_line.alt[i] += allele_addition
				end
			end
			vcf_line.sample_names.each do |sample_name|
				sample1 = prev_line.samples[sample_name]
				sample2 = vcf_line.samples[sample_name]

				geno1 = sample1[:"GT"]
				geno1 = './.' if geno1.nil?
				geno1_ar = geno1.split(/[\/\|\\]/)
				geno2 = sample2[:"GT"]
				geno2 = './.' if geno2.nil?
				geno2_ar = geno2.split(/[\/\|\\]/)
				# adjust allele index 
				geno2_ar.map! do |allele_i| 
					if allele_i == '0' || allele_i == '.'
						allele_i
					else
						allele_i.to_i + prev_line.alt.length  
					end
				end
				geno2 = geno2_ar.join("/")
				new_geno = Array.new

				# get genotype qualities
				gq1 = sample1[:"GQ"]
				gq2 = sample2[:"GQ"]
				if(gq1.nil?)
					gq1=-1
				else
					gq1=gq1.to_i
				end
				if(gq2.nil?)
					gq2=-1
				else
					gq2=gq2.to_i
				end
				
				if(geno1 =~ /[1-9]/ && geno2 =~ /[1-9]/) # both are variant
					STDERR.print "WARNING: multiple variants were called in the same sample (#{sample_name}) at the same site #{vcf_line.chr}:#{vcf_line.pos} GQs: #{gq1} #{gq2} "
					#if(geno1 =~ /[0\.]/ && geno2 =~ /[0\.]/) # both are heterozygous, so they can be merged
					#	if(geno1 =~ /[1-9]/)
					#		new_geno.push $&[0]
					#	end
					#	if(geno2 =~ /[1-9]/)
					#		new_geno.push $&[0]
					#	end
					#	STDERR.puts "both are heterozygous, merging to #{new_geno}"
					#	prev_line.samples[sample_name][:"GT"] = new_geno.join("/")
					if(gq1 >= gq2)
						STDERR.puts "keeping the first genotype"
						# nothing needs to be done
					else
						STDERR.puts "keeping the second genotype"
						prev_line.samples[sample_name] = vcf_line.samples[sample_name]
						prev_line.samples[sample_name][:"GT"] = geno2
					end
				elsif(geno2 =~ /[1-9]/ || geno1 == './.') # geno2 is variant or geno1 is nil
					prev_line.samples[sample_name] = vcf_line.samples[sample_name]
					prev_line.samples[sample_name][:"GT"] = geno2
				elsif(geno1_ar[0]=='0' && geno1_ar[1] =='0' && geno2_ar[0] == '0' && geno2_ar[1]=='0') # all are ref
					begin
						if(sample1[:VR].to_i >= sample2[:VR].to_i)
							# nothing needs to be done
						else
							prev_line.samples[sample_name] = vcf_line.samples[sample_name]
						end
					rescue
						STDERR.puts "No variant read data available, will just keep first genotype when both ref"
					end
				elsif(geno1 =~ /[1-9]/ || geno2 == './.')
					# nothing needs to be done
				else
					raise "Unhandled case: geno1 #{geno1}    geno2 #{geno2}"
				end


					
			end
			prev_line.alt += vcf_line.alt
		else
			puts prev_line
			prev_line = vcf_line
		end
	else
		prev_line = vcf_line
	end
end
puts prev_line.print() unless prev_line.nil?
STDERR.puts "FINISHED, merged #{merge_count} rows"
