#!/usr/bin/ruby
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 1)
	puts "USAGE: ruby split_multiallelic_rows.rb [vcf_file] [-info]"
	exit 0
end

labels=nil
prev_line = nil	
split_count = 0
if("-info" == ARGV[1])
	@info=true
end
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
#	if(vcf_line.alt.length<2)
#		puts vcf_line.print()
#		next
#	end
	vcf_line.alt().each_with_index do |alt, alt_i|
		new_vcf_line = Vcf_line.read_line(line, false, labels)
		new_ref, new_alt = Vcf_line.simplify(new_vcf_line.ref.to_s,alt.to_s)
		new_vcf_line.ref=new_ref
		new_vcf_line.alt=[new_alt]
		if(@info)
			new_vcf_line.info.each_key do |info_key|
				if(new_vcf_line.info[info_key].class == Array)
					new_vcf_line.info[info_key] = [vcf_line.info[info_key][alt_i]]
				end
			end
		end
		new_vcf_line.sample_names().each do |sample_name|
			geno = new_vcf_line.samples()[sample_name][:GT].to_s
			if(geno =~ /[\/|]/)
				sep = $&
			else
				sep = '/'
			end
			sample_alleles = geno.split(/[\/|]/)
			0.upto(sample_alleles.length).each do |sample_allele_i|
				if(sample_alleles[sample_allele_i].to_i == alt_i+1)
					sample_alleles[sample_allele_i] = 1
				elsif(sample_alleles[sample_allele_i] =~ /[1-9]/)
					sample_alleles[sample_allele_i] = '.' # sample is some other variant allele
				end
			end
			if(sample_alleles.length == 0)
				new_vcf_line.samples[sample_name][:GT] = './.'
			else	
				new_vcf_line.samples[sample_name][:GT] = sample_alleles.join(sep)
			end
		end
		puts new_vcf_line.print()
	end
	split_count+=1 if vcf_line.alt.length > 1
end

STDERR.puts "FINISHED, split #{split_count} rows"
