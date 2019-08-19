#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Ss_metric_extractor

	def initialize()
		@sample_metrics = Array.new
	end

	def run
		labels=nil
		STDIN.each_line do |line|
			line.chomp!
			if(line[0,1] == '#')
				if( line[0,4] == '#CHR' )
					labels = line.split("\t") 
				end
				next
			end
			vcf_line = Vcf_line.read_line(line, false, labels)
			process_line(vcf_line)
		end
		print_results
	end

	def process_line(line)
		line.sample_names.each_with_index do |sample_name, sample_index|
			return if line.filter.to_s != "PASS"

			#set up for sample
			if(@sample_metrics[sample_index].nil?)
				@sample_metrics[sample_index] = Hash.new(0) 
				@sample_metrics[sample_index][:name] = sample_name
				@sample_metrics[sample_index][:lengths] = Hash.new(0) # store a count of indel lengths for each sample
				@sample_metrics[sample_index][:AF] = Hash.new(0) # store a count of allele frequency bins
			end

			# frameshift_count
			if line.samples[sample_name][:"GT"].to_s =~ /([0-9.])\/([0-9.])/
				a1=$1.to_i
				a2=$2.to_i
			else
				next
			end
			next if(a1 == 0 && a2 == 0)
			# indel_count
			@sample_metrics[sample_index][:indel_count] += 1 
			allele_codes = [a1,a2]
			allele_codes.delete(0)
			if(allele_codes.length == 1)
				@sample_metrics[sample_index][:heter_count] += 1
				heter=true
			elsif(allele_codes.length == 2 && allele_codes[0] != allele_codes[1])
				@sample_metrics[sample_index][:heter_count] += 2
				heter=true
			else
				heter=false
			end
			allele_codes.uniq!
			allele_codes.each do |allele_code|
				allele = line.alt[allele_code-1].to_s
				raise "bad allele!\n#{line.print}" if allele == ""
				length = allele.length - line.ref.to_s.length
				raise "not an indel!\n#{line.print}" if length == 0
				# frameshift count
				if(length.abs % 3 != 0)
					@sample_metrics[sample_index][:frameshift_count] += 1 
				end
				# insertion_count
				if(length > 0)
					@sample_metrics[sample_index][:insertion_count] += 1 
				end
				# 1bp indel count
				if(length.abs == 1)
					@sample_metrics[sample_index][:count_1bp] += 1 
				end
				# greater than 10bp and indel length counts
				if(length.abs > 10)
					@sample_metrics[sample_index][:count_gt10] += 1 
					if(length > 0)
						@sample_metrics[sample_index][:lengths][:">10"] += 1 
					else
						@sample_metrics[sample_index][:lengths][:"<-10"] += 1 
					end
				else
					@sample_metrics[sample_index][:lengths][length] += 1 
				end
			#	if(line.info[:"Mask"] < 6)
					@sample_metrics[sample_index][:count_low_complexity] += 1
			#	end
				if(line.info[:MAXIMPACT] == :"HIGH")
					@sample_metrics[sample_index][:count_high_impact] += 1
				end
				if(line.info[:MAXIMPACT] == :"MODERATE")
					@sample_metrics[sample_index][:count_moderate_impact] += 1
				end
				if(line.info[:"pcrSlip"] == :"T")
					@sample_metrics[sample_index][:count_slippage] += 1
				end
				if(line.info[:"repeatLength"] == 1)
					@sample_metrics[sample_index][:count_homo] += 1
				end
				if(line.info[:"repeatLength"] != nil && line.info[:"repeatLength"] != :"." && line.info[:"repeatLength"] > 1)
					@sample_metrics[sample_index][:count_repeats] += 1
				end
				aaf = line.info[:AAF]
				ac = line.info[:AC]
				if(!ac.nil? && ac.to_i==1)
					@sample_metrics[sample_index][:singleton] += 1
				end
				unless(aaf.nil?)
					aaf = aaf.to_f
					maf = aaf
					maf = 1.0-aaf if aaf > 0.5
					if(aaf < 0.0003162278)
						@sample_metrics[sample_index][:AF][:"af0.0003"] += 1
					elsif(aaf < 0.001)
						@sample_metrics[sample_index][:AF][:"af0.001"] += 1
					elsif(aaf < 0.0031622779)
						@sample_metrics[sample_index][:AF][:"af0.003"] += 1
					elsif(aaf < 0.01)
						@sample_metrics[sample_index][:AF][:"af0.01"] += 1
					elsif(aaf < 0.0316227781)
						@sample_metrics[sample_index][:AF][:"af0.03"] += 1
					elsif(aaf < 0.1)
						@sample_metrics[sample_index][:AF][:"af0.1"] += 1
					elsif(aaf < 0.3162277709)
						@sample_metrics[sample_index][:AF][:"af0.3"] += 1
					else
						@sample_metrics[sample_index][:AF][:"af1.0"] += 1
					end
						
				end
				#	"	#{sample[:AF][:"af0.0003"]}	#{sample[:AF][:"af0.001"]}	#{sample[:AF][:"af0.003"]}	#{sample[:AF][:"af0.01"]}	#{sample[:AF][:"af0.03"]}	#{sample[:AF][:"af0.1"]}	#{sample[:AF][:"af0.3"]}	#{sample[:AF][:"af1.0"]}"
				if(line.info[:effect].to_s.include?("damaging") || line.info[:provean_prediction].to_s.include?("Deleterious"))
					@sample_metrics[sample_index][:count_damaging] += 1
				end

				#	unless(heter)
				#		@sample_metrics[sample_index][:count_homoz_damaging] += 1
				#	end
				#	
				#end
			end
		end
	end

	def print_results
		puts "#sample	indel_count	frameshift_count	frameshift_rate	insertion_count	deletion_count	del_ins_ratio	count_1bp	ratio_1bp	gt10bp_count	gt10bp_ratio	heter_count	heter_ratio	low_complexity_count	low_complexity_ratio	high_impact_count	high_impact_ratio	moderate_impact_count	moderate_impact_ratio	slippage_count	slippage_ratio	homopoly_count	homopoly_ratio	repeat_count	repeat_ratio	damaging_count	damaging_ratio	af0.0003	af0.001	af0.003	af0.01	af0.03	af0.1	af0.3	af1.0	singleton	<-10	#{(-10..10).to_a.join("\t")}	>10"
		@sample_metrics.each do |sample|
			deletion_count = sample[:indel_count]-sample[:insertion_count]
			print "#{sample[:name]}	#{sample[:indel_count]}	#{sample[:frameshift_count]}	#{sample[:frameshift_count].to_f/sample[:indel_count].to_f}	#{sample[:insertion_count]}	#{deletion_count}	#{deletion_count.to_f/sample[:insertion_count].to_f}	#{sample[:count_1bp]}	#{sample[:count_1bp].to_f/sample[:indel_count].to_f}	#{sample[:count_gt10]}	#{sample[:count_gt10].to_f/sample[:indel_count].to_f}	#{sample[:heter_count]}	#{sample[:heter_count].to_f/sample[:indel_count].to_f}	#{sample[:count_low_complexity]}	#{sample[:count_low_complexity].to_f/sample[:indel_count].to_f}	#{sample[:count_high_impact]}	#{sample[:count_high_impact].to_f/sample[:indel_count].to_f}	#{sample[:count_moderate_impact]}	#{sample[:count_moderate_impact].to_f/sample[:indel_count].to_f}	#{sample[:count_slippage]}	#{sample[:count_slippage].to_f/sample[:indel_count].to_f}	#{sample[:count_homo]}	#{sample[:count_homo].to_f/sample[:indel_count].to_f}	#{sample[:count_repeats]}	#{sample[:count_repeats].to_f/sample[:indel_count].to_f}	#{sample[:count_damaging]}	#{sample[:count_damaging].to_f/sample[:indel_count].to_f}	#{sample[:AF][:"af0.0003"]}	#{sample[:AF][:"af0.001"]}	#{sample[:AF][:"af0.003"]}	#{sample[:AF][:"af0.01"]}	#{sample[:AF][:"af0.03"]}	#{sample[:AF][:"af0.1"]}	#{sample[:AF][:"af0.3"]}	#{sample[:AF][:"af1.0"]}	#{sample[:singleton]}"
			print "\t#{sample[:lengths][:"<-10"]}"
			(-10..10).to_a.each do |length|
				print "\t#{sample[:lengths][length]}"
			end
			puts "\t#{sample[:lengths][:">10"]}"
		end
	end
end

if(ARGV[0] == "-h" || ARGV[0] == "--help")
	puts "USAGE: cat vcf_file | ruby extract_single_sample_metrics.rb"
	exit 1
end

extractor = Ss_metric_extractor.new()
extractor.run
