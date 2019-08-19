#!/usr/bin/ruby
#merge 2 VCF files giving the first priority with all conflicts (however it will always keep the GL if one of them has a nil GL)
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 2)
	puts "USAGE: ruby copy_sample_entry.rb [main.vcf.gz] [source.vcf.gz] comma,sep,list,of,entries"
	exit 0
end

vcf_file1=ARGV[0]
vcf_file2=ARGV[1]
vcf_reader1=IO.popen("zgrep -v '^#' #{vcf_file1}")
vcf_reader2=IO.popen("zgrep -v '^#' #{vcf_file2}")
entries=ARGV[2].split(",")
entries.map! {|key| key.to_sym}
source1 = "1" if source1.nil?
source2 = "2" if source2.nil?


labels1= `zgrep -m1 "#CHROM	" #{vcf_file1}`.chomp.split("\t")
labels2= `zgrep -m1 "#CHROM	" #{vcf_file2}`.chomp.split("\t")

#`zcat #{vcf_file2} |cut -f 1,2  | grep -v '#' > tmp.#{$$}.sites`

# print header
puts `zcat #{vcf_file1} | head -n 300 | grep '^#' `.chomp

#File.open("tmp.#{$$}.sites",'r').each_line do |pos_str|
begin
	while(true)
		line1,line2,vcf1,vcf2,merged_vcf,source = nil
		line1 = vcf_reader1.readline.chomp
		line2 = vcf_reader2.readline.chomp
		# find and load VCF lines
		#line1=`tabix #{vcf_file1} #{chr}:#{coor}-#{coor} | grep "	#{coor}	"`.chomp
		#line2=`tabix #{vcf_file2} #{chr}:#{coor}-#{coor} | grep "	#{coor}	"`.chomp
		vcf1 = Vcf_line.read_line(line1, false, labels1)
		vcf2 = Vcf_line.read_line(line2, false, labels2) 
		entries.each do |key|
			vcf1.format.push(key) unless vcf1.format.include?(key)
		end
		vcf1.sample_names.each do |sample_name|
			entries.each do |key|
				vcf1.samples[sample_name][key]=vcf2.samples[sample_name][key]
			end
		end
		puts vcf1.print
		STDOUT.flush
	end
rescue EOFError => err
	vcf_reader1.close
	vcf_reader2.close
end


`rm -f tmp.#{$$}.sites`
STDERR.puts "FINISHED"
