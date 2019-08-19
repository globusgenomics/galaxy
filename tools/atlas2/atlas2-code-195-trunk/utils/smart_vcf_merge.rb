#!/usr/bin/ruby
#merge 2 VCF files giving the first priority with all conflicts (however it will always keep the GL if one of them has a nil GL)
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 2)
	puts "USAGE: ruby smart_vcf_merge.rb [priority_vcf.gz] [vcf_file2.gz] [source1] [source2]"
	exit 0
end

vcf_file1=ARGV[0]
vcf_file2=ARGV[1]
source1=ARGV[2]
source2=ARGV[3]
source1 = "1" if source1.nil?
source2 = "2" if source2.nil?


labels1= `zgrep -m1 "#CHROM	" #{vcf_file1}`.chomp.split("\t")
labels2= `zgrep -m1 "#CHROM	" #{vcf_file2}`.chomp.split("\t")
sample_names = labels1[9..(labels1.length)] + labels2[9..(labels2.length)]
sample_names.sort!
sample_names.uniq!

`zcat #{vcf_file1} #{vcf_file2} |cut -f 1,2  | grep -v '#' | sort | uniq > tmp.#{$$}`
`grep -P '^[0-9]' tmp.#{$$} | sort -k1,1n -k2,2n > tmp.#{$$}.sites`
`grep -vP '^[0-9]' tmp.#{$$} | sort -k1,1d -k2,2n >> tmp.#{$$}.sites`
`rm -f tmp.#{$$}`

# print header
puts `zcat #{vcf_file1} | head -n 300 | grep '##' `.chomp
puts "#{labels1[0..8].join("\t")}	#{sample_names.join("\t")}"

File.open("tmp.#{$$}.sites",'r').each_line do |pos_str|
	# find and load VCF lines
	chr,coor,line1,line2,vcf1,vcf2,merged_vcf,source = nil
	chr,coor = pos_str.chomp.split("\t")	
	line1=`tabix #{vcf_file1} #{chr}:#{coor}-#{coor} | grep "	#{coor}	"`.chomp
	line2=`tabix #{vcf_file2} #{chr}:#{coor}-#{coor} | grep "	#{coor}	"`.chomp
	if(line1.lines.count > 1 || line2.lines.count >1)
		raise "One of the VCF files has multiple lines with the same coordinate, that's not allowed by this program"
	end

	if(line1!='' && line2!='') # both file have a VCF line here
		vcf1 = Vcf_line.read_line(line1, false, labels1)
		vcf2 = Vcf_line.read_line(line2, false, labels2) 
		merged_vcf = Vcf_line.read_line(line1, false, labels1)
		merged_vcf.sample_names = sample_names
		merged_vcf.samples = Hash.new
		source=[source1,source2]
		

		# merge alt column
		allele_hash1=Hash.new # keep track of original allele indices
		allele_hash2=Hash.new
		allele_addition1 = ''
		allele_addition2 = ''
		# select longest ref and added needed bases to ALTs
		if(vcf1.ref.to_s.length > vcf2.ref.to_s.length)
			allele_addition2 = vcf1.ref.to_s[(vcf2.ref.to_s.length)...(vcf1.ref.to_s.length)]
			merged_vcf.ref = vcf1.ref
		elsif(vcf2.ref.to_s.length > vcf1.ref.to_s.length)
			allele_addition1 =  vcf2.ref.to_s[(vcf1.ref.to_s.length)...(vcf2.ref.to_s.length)]
			merged_vcf.ref = vcf2.ref
		end
		vcf1.alt.each_with_index do |allele, i|
			vcf1.alt[i] = allele.to_s
			vcf1.alt[i] += allele_addition1
		end
		vcf2.alt.each_with_index do |allele, i|
			vcf2.alt[i] = allele.to_s
			vcf2.alt[i] += allele_addition2
		end
		# do actual merging of alts
		merged_vcf.alt = (vcf1.alt + vcf2.alt).sort.uniq

		# track original and new allele index
		vcf1.alt.each_with_index do |allele, i|
			allele_hash1[i+1] = merged_vcf.alt.index(allele)+1
		end
		vcf2.alt.each_with_index do |allele, i|
			allele_hash2[i+1] = merged_vcf.alt.index(allele)+1
		end

		# merge sample columns
		sample_names.each do |name|
			sample1,sample2,geno1,gl1,gl2,geno,new_index=nil
			sample1 = vcf1.samples[name]
			sample2 = vcf2.samples[name]
			geno1 = sample1[:GT] unless sample1.nil?
			gl1 = sample1[:GL] unless sample1.nil?
			gl2 = sample2[:GL] unless sample2.nil?
			if(geno1.nil? || geno1 == './.' || geno1 == '.')
				if( sample2.nil? || sample2[:GT].nil?)  # both files have no geno
					merged_vcf.samples[name]=sample1
					merged_vcf.samples[name][:GL] = gl2 if gl1.nil?
				else # replace geno index
					geno = sample2[:GT].to_s
					if(geno =~ /[1-9]+/)
						new_index=allele_hash2[$&.to_i].to_s
						vcf2.samples[name][:GT]=geno.gsub($&[0], new_index)
					end
					merged_vcf.samples[name]=sample2
					merged_vcf.samples[name][:GL] = gl1 if gl2.nil?
				end
			else
				geno = sample1[:GT].to_s
				if(geno =~ /[1-9]+/)
					new_index=allele_hash1[$&.to_i].to_s
					vcf1.samples[name][:GT]=geno.gsub($&[0], new_index)
				end
				merged_vcf.samples[name]=sample1
				merged_vcf.samples[name][:GL] = gl2 if gl1.nil?
			end
		end
	else # only one VCF file has a line here
		if(line1.nil? || line1=="")
			merged_vcf = Vcf_line.read_line(line2, false, labels2)
			source=source2
		else
			merged_vcf = Vcf_line.read_line(line1, false, labels1)
			source=source1
		end
		merged_vcf.sample_names = sample_names
	end

	merged_vcf.format=[:GT,:GL]
	merged_vcf.info=Hash.new()	
	merged_vcf.info[:SF]=source
	puts merged_vcf.print
end


`rm -f tmp.#{$$}.sites`
STDERR.puts "FINISHED"
