#!/usr/bin/ruby


$:.unshift File.join(File.dirname(__FILE__),'.')
require 'ref_seq.rb'


if(ARGV.length < 4 )
	puts "USAGE: ruby generate_random_indels.rb ref.fasta target.bed length_list #sites"
	exit 1
end


ref_seq = Ref_seq.new(ARGV[0])
num_sites = ARGV[3].to_i

site_list = Array.new # list of genome positions that will be randomly sampled from
length_list = Array.new # list of indel lengths that will be randomly sampled from
bases = ['A','T','G','C']

# read length list
STDERR.puts "reading indel length list"
File.open(ARGV[2],'r').each_line do |length|
	length.chomp!
	length_list.push(length.to_i)
end


# read in BED regions
STDERR.puts "creating loci pool from BED file"
File.open(ARGV[1], 'r').each_line do |entry|
	cols = entry.split("\t")
	chr = cols[0]
	start = cols[1].to_i + 1
	stop = cols[2].to_i
	start.upto(stop) do |coor|
		site_list.push([chr,coor])
	end
end

puts "##fileformat=VCFv4.1"
puts "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
count = 0

STDERR.puts "Generating random indels"
while(count < num_sites)
	chr,coor = site_list[rand(site_list.length)]
	length = length_list[rand(length_list.length)]
	
	ref = ""
	alt = ""
	if(length>0) # insersion
		ref = ref_seq.get_base(chr, coor)
		alt = ref
		while(alt.length-1 < length)
			alt += bases[rand(bases.length)] # insert random bases
		end
	elsif(length<0)
		ref = ref_seq.get_base_range(chr, coor, coor+length.abs)
		alt = ref_seq.get_base(chr, coor)
	else
		next # this is not an indel length	
	end
	count += 1
	puts "#{chr}	#{coor}	.	#{ref}	#{alt}	.	.	."
end
