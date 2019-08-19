#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'ref_seq.rb'
require 'vcf_line.rb'





def getline
	begin
		line = $vcf_reader.readline
		line.chomp!
	rescue EOFError
		line = :eof
	end
	return line
end





if(ARGV.length < 4)
	puts "USAGE: update_ref_seq.rb <original_ref.fasta> <SNPs.vcf> <output.fasta> <changes.vcf>"
	exit 1	
end

outfile = ARGV[2]
if(File.exists?(outfile))
	STDERR.puts "The output file #{outfile} already exists, please rename or choose a different file"
	exit 1
end
change_writer = File.open(ARGV[3], 'w')

print "Making chromosome list..."
fastq_headers = []
chromosome_list = []
open(ARGV[0], 'r').grep(/^>([a-zA-Z0-9_\.\-]+).*/) do |chr|
	fastq_headers.push(chr)
	chromosome_list.push($1.to_sym)
end
puts "DONE"

ref = Ref_seq.new(ARGV[0])

$vcf_reader = File.open(ARGV[1], 'r')
line = getline()
labels = ""
while(line[0]=='#')
	if(line[0,6] == '#CHROM')
		labels = line.split("\t")
		change_writer.print line
	end
	line = getline()
end
if(line == :eof)
	STDERR.puts "The VCF file has no data!"
	exit 1
end
change_count = 0
vcf_line = Vcf_line.read_line(line, false, labels) 
chromosome_list.each_with_index do |current_chr, chr_i|
	if(current_chr != ref.current_chr)
		ref.load_next_chr(current_chr)
	end
	#if(line == :eof || chromosome_list.index(vcf_line.chr.to_sym) > chr_i)
	#	print "Writing new chromosome #{current_chr}..."
	#	ref.write_chromosome_to_fasta(outfile, fastq_headers[chr_i])
	#	puts "DONE"
	#	next
	if(line != :eof)
		if(chromosome_list.index(vcf_line.chr.to_sym) < chr_i)
			STDERR.puts "ERROR: The reference genome and the VCF file do not have their chromosomes in the same order."
			exit 1
		end
		while(vcf_line.chr.to_sym == current_chr)
			max_allele_i = nil
			max_af = 0
			sum_af = 0.0
			if(vcf_line.info[:"AC"].nil? || vcf_line.info[:"AN"].nil?)
				STDERR.puts "ERROR: The VCF file does not include AN and/or AC in the INFO at line #{current_chr}:#{vcf_line.pos}"
				exit 1
			end
			an = vcf_line.info[:"AN"]
			vcf_line.alt.each_with_index do |allele, allele_i|
				if(vcf_line.info[:"AC"].class == Array)
					ac = vcf_line.info[:"AC"][allele_i]
				else
					ac= vcf_line.info[:"AC"]
				end
				af = ac.to_f / an.to_f # calc allele frequency
				sum_af += af
				next if(vcf_line.ref.length != allele.length)
				if( af > max_af )
					max_af = af
					max_allele_i = allele_i
				end
			end
			# update the reference if the highest alt AF > the REF AF
			if(max_af > 1-sum_af && an > 300)
				raise "alt allele seems to be same as ref: #{vcf_line.ref}->#{vcf_line.alt[max_allele_i]}" if vcf_line.ref[0] == vcf_line.alt[max_allele_i][0] # sanity check
				ref.change_base(current_chr, vcf_line.pos, vcf_line.alt[max_allele_i][0])
				change_count += 1
				change_writer.puts vcf_line.print
				print "Changes: #{change_count}            \r"
			end

			line = getline()
			if(line == :eof)
				break
			else
				vcf_line = Vcf_line.read_line(line, false, labels)
			end
		end
	end
	print "Writing new chromosome #{current_chr}..."
	ref.write_chromosome_to_fasta(outfile, fastq_headers[chr_i])
	puts "DONE"
end


$vcf_reader.close
change_writer.close
puts "FINISHED"
puts "Changed the reference at #{change_count} sites"

