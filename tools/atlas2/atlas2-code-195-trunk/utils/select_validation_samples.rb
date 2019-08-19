$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 2)
	puts "USAGE: ruby select_validation_samples.vcf source_vcf sites sample_inventory #samples"
	exit 0
end

source_vcf = ARGV[0]
labels=nil
source_labels = `grep -m1 "#CHROM" #{source_vcf}`.chomp.split("\t")
sample_inventory = Array.new
num_samples = ARGV[3].to_i
File.open(ARGV[2], 'r').each_line do |line|
	sample_inventory.push(line.chomp)
end
File.open(ARGV[1], 'r').each_line do |line|
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
			puts source_labels.join("\t")
		end
		next
	end
	picked_samples = Array.new
	sites_vcf_line = Vcf_line.read_line(line, false, labels)
	source_line = `grep -m 1 "#{sites_vcf_line.chr}	#{sites_vcf_line.pos}	" #{source_vcf}`.chomp
	raise "could not find site #{sites_vcf_line.chr}:#{sites_vcf_line.pos} in the source_vcf" if source_line.nil? || source_line == ''
	vcf_line = Vcf_line.read_line(source_line, false, source_labels)
	variant_samples = Array.new
	vcf_line.samples.each do |name, data|
		if(data[:"GT"].to_s =~ /[1-9]/ && sample_inventory.include?(name))
			variant_samples.push(name)
		end
	end
	while(picked_samples.length < num_samples && variant_samples.length > 0)
		#if(variant_samples.length <= num_samples)
		#	picked_samples += variant_samples
		#	variant_samples = Array.new
		#end
		if(variant_samples.length > 0)
			i = rand(variant_samples.length)
			picked_samples.push(variant_samples.delete_at(i)) 
		#else
		#	sample_name = vcf_line.sample_names[rand(vcf_line.sample_names.length)]
		#	if( (!picked_samples.include?(sample_name)) && sample_inventory.include?(sample_name))
		#		picked_samples.push(sample_name)
		#	end
		end
	end
	vcf_line.info[:"samples"]=picked_samples
	#vcf_line.samples=Hash.new
	#vcf_line.sample_names=Array.new
	puts vcf_line.print
end
