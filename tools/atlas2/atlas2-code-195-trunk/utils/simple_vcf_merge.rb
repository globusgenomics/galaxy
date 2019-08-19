#!/usr/bin/ruby
#merge 2 VCF files 
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'


def fill_blank_samples(vcf, sample_names)
	vcf.sample_names=sample_names
	sample_names.each do |name|
		if(vcf.samples[name].nil?)
			vcf.samples[name]=Hash.new
			vcf.samples[name][:GT]='./.'
			vcf.samples[name][:GL]='.,.,.'
		elsif(vcf.samples[name][:GL]=='.')
			vcf.samples[name][:GL]='.,.,.'
		end
	end
end


if(ARGV.length < 2)
	puts "USAGE: ruby simple_vcf_merge.rb [vcf_file1.gz] [vcf_file2.gz] [source1] [source2]"
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
puts `zcat #{vcf_file1} | head -n 200 | grep '##' `.chomp
puts '##INFO=<ID=VT,Number=1,Type=String,Description="indicates what type of variant the line represents">'
puts "#{labels1[0..8].join("\t")}	#{sample_names.join("\t")}"

File.open("tmp.#{$$}.sites",'r').each_line do |pos_str|
	# find and load VCF lines
	chr,coor,lines1,lines2,vcf1,vcf2= nil
	chr,coor = pos_str.chomp.split("\t")	
	lines1=`tabix #{vcf_file1} #{chr}:#{coor}-#{coor} | grep "	#{coor}	"`.chomp.split("\n")
	lines2=`tabix #{vcf_file2} #{chr}:#{coor}-#{coor} | grep "	#{coor}	"`.chomp.split("\n")
	
	lines1.each do |vcf_str|
		vcf1 = Vcf_line.read_line(vcf_str, false, labels1)
		fill_blank_samples(vcf1, sample_names)
		vcf1.format=[:GT,:GL]
		vcf1.info[:VT]=source1
		puts vcf1.print
	end
		
	lines2.each do |vcf_str|
		vcf2 = Vcf_line.read_line(vcf_str, false, labels2)
		fill_blank_samples(vcf2, sample_names)
		vcf2.format=[:GT,:GL]
		vcf2.info[:VT]=source2
		puts vcf2.print
	end

end


`rm -f tmp.#{$$}.sites`
STDERR.puts "FINISHED"
