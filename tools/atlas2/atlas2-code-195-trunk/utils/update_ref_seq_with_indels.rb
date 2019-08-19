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



def write_chromosome_to_fasta( filename, header_line=">#{@chromosome}", bases="" )
	writer = File.open(filename, 'a')
	writer.puts ">#{header_line}";
	(0...(bases.length)).step(60) do |i|
		writer.puts bases[i...(i+60)].join('')
	end
	writer.close
end





if(ARGV.length < 5)
	puts "USAGE: update_ref_seq.rb <original_ref.fasta> <variants.vcf> <output.fasta> <output.chain> <change_list.vcf>"
	exit 1	
end

outfile = ARGV[2]
chain_writer = File.open(ARGV[3], 'w')
rev_chain_writer = File.open("#{ARGV[3]}.reverse", 'w')
change_writer = File.open(ARGV[4], 'w')
if(File.exists?(outfile))
	STDERR.puts "The output file #{outfile} already exists, please rename or choose a different file"
	exit 1
end

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
	chain_rows = Array.new
	chain_size = 0
#	ref.load_next_chr(current_chr)
	old_ref = `samtools faidx #{ARGV[0]} #{current_chr}`.sub(/>[chr0-9]*/,'').gsub(/\n/,'')
	old_ref_i = 0
	new_ref = Array.new
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
		while(!vcf_line.nil? && vcf_line.chr.to_sym == current_chr)
			while(vcf_line.pos > old_ref_i + 1)
				new_ref.push(old_ref[old_ref_i])
				old_ref_i += 1
				chain_size += 1
			end
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
				# next if(vcf_line.ref.length != allele.length)
				if( af > max_af )
					max_af = af
					max_allele_i = allele_i
				end
			end
			# update the reference if the highest alt AF > the REF AF (also skip large deletions) - exclude indel with length > 3000
			if(max_af > 1-sum_af && vcf_line.alt[max_allele_i] != :"<DEL>")
				raise "alt allele seems to be same as ref: #{vcf_line.ref}->#{vcf_line.alt[max_allele_i]}" if vcf_line.ref == vcf_line.alt[max_allele_i] # sanity check
				if(vcf_line.ref.length == vcf_line.alt[max_allele_i].length) # SNP
					new_ref.push(vcf_line.alt[max_allele_i][0])
					old_ref_i += 1
					chain_size += 1
				elsif(vcf_line.ref.length > vcf_line.alt[max_allele_i].length && vcf_line.ref.length - vcf_line.alt[max_allele_i].length < 3000) # deletion
					new_ref.push(old_ref[vcf_line.pos-1]) # the first base is still ref
					indel_length = vcf_line.ref.length - vcf_line.alt[max_allele_i].length
					old_ref_i += (indel_length + 1)
					chain_size += 1
					chain_rows.push [chain_size,indel_length,0]
					chain_size = 0
				elsif( vcf_line.alt[max_allele_i].length - vcf_line.ref.length < 3000) # insertion
					new_ref.push(old_ref[vcf_line.pos-1]) # the first base is still ref
					old_ref_i += 1
					chain_size += 1
					ins = vcf_line.alt[max_allele_i].to_s.sub(/#{vcf_line.ref}/, "")
					ins.each_char do |base|
						new_ref.push(base)
					end
					chain_rows.push [chain_size,0,ins.length]
					chain_size = 0
				else
					raise "unhandled case"
				end
				change_count += 1
				print "Changes: #{change_count}            \r"
				change_writer.puts vcf_line.print
			else # keep reference
				new_ref.push(old_ref[vcf_line.pos-1])
				old_ref_i += 1
				chain_size += 1
			end

			while(!vcf_line.nil? && vcf_line.chr.to_sym == current_chr && vcf_line.pos <= old_ref_i + 1) # get next line, may have to skip a few if deletion spans them
				line = getline()
				if(line == :eof)
					vcf_line = nil
					break
				else
					vcf_line = Vcf_line.read_line(line, false, labels)
				end
			end
		end # while(vcf_line.chr.to_sym == current_chr)
	end # if(line != :eof)
	while(old_ref_i < old_ref.length) # add the end of the chromosome
		new_ref.push(old_ref[old_ref_i])
		old_ref_i += 1
		chain_size += 1
	end
	print "Writing new chromosome #{current_chr}..."
	write_chromosome_to_fasta(outfile, current_chr, new_ref)
	puts "DONE"
	print "Writing chromosome #{current_chr} chain file entry..."
	chain_writer.puts "chain 100 #{current_chr} #{old_ref.length} + 1 #{old_ref.length} #{current_chr} #{new_ref.length} + 1 #{new_ref.length} #{chr_i+1}"
	rev_chain_writer.puts "chain 100 #{current_chr} #{new_ref.length} + 1 #{new_ref.length} #{current_chr} #{old_ref.length} + 1 #{old_ref.length} #{chr_i+1}"
	chain_rows.each do |row|
		chain_writer.puts row.join("	")
		rev_chain_writer.puts "#{row[0]}	#{row[2]}	#{row[1]}"
	end
	chain_writer.puts chain_size - 1
	chain_writer.puts ""
	rev_chain_writer.puts chain_size - 1
	rev_chain_writer.puts ""
	puts "DONE"
end # for each chromosome


$vcf_reader.close
chain_writer.close
rev_chain_writer.close
change_writer.close
puts "FINISHED"
puts "Changed the reference at #{change_count} sites"

