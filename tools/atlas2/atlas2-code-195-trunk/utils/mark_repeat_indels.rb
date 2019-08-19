#!/usr/bin/ruby

def get_repeat_substr(seq)
	if(seq =~ /(.+?)\1+/) # find shortest repeat unit
		if($&.length == seq.length)
			return $1
		end
	end
	return seq
end

if __FILE__ == $0
	$:.unshift File.join(File.dirname(__FILE__),'.')
	require 'vcf_line'
	require 'ref_seq'
	window_size = 100



	if(ARGV.length < 2)
		puts "USAGE: ruby mark_repeat_indels.rb vcf_file reference_genome"
		exit(0)
	end

	vcf_file = ARGV[0]
	ref_genome = ARGV[1]
	reference = Ref_seq.new(ref_genome)
	labels=nil
	File.open(vcf_file, 'r').each_line do |line|
		if(line[0,1] == '#')
			if( line[0,4] == '#CHR' )
				labels = line.split("\t") 
				puts "##INFO=<ID=repeat,Number=A,Type=Character,Description=\"Is the allele a repeat indel? (T/F)\">"
				puts "##INFO=<ID=repeatStr,Number=A,Type=String,Description=\"If allele is a repeat indel, this indicates the duplicated or deleted repeat sequence\">"
				puts "##INFO=<ID=repeatCnt,Number=A,Type=Integer,Description=\"If allele is a repeat indel, this indicates the number of times the repeatStr is repeated\">"
				puts "##INFO=<ID=repeatLength,Number=A,Type=Integer,Description=\"If allele is a repeat indel, this indicates the sequence length of the repeatStr\">"
				puts "##INFO=<ID=repeatRegLength,Number=A,Type=Integer,Description=\"If allele is a repeat indel, this indicates the sequence length of the repeat region\">"
			end
			puts line
			next
		end
		vcf_line = Vcf_line.read_line(line, false, labels)
		#next if vcf_line.alt.length > 1 # skip multi-allelic sites
		vcf_line.info[:"repeat"]=Array.new
		vcf_line.info[:"repeatStr"]=Array.new
		vcf_line.info[:"repeatCnt"]=Array.new
		vcf_line.info[:"repeatLength"]=Array.new
		vcf_line.info[:"repeatRegLength"]=Array.new
		vcf_line.alt.each do |allele|
			ref = vcf_line.ref.to_s
			alt= allele.to_s
			ref,alt = Vcf_line.simplify(vcf_line.ref.to_s,allele.to_s) if ref.length != 1 && alt.length != 1
			if(alt.length < ref.length) # deletion
				indel_seq = ref.dup
				indel_seq[0]=''
				rep_str = get_repeat_substr(indel_seq)
				leading = reference.get_base_range(vcf_line.chr, vcf_line.pos-rep_str.length+1, vcf_line.pos)
				trailing = reference.get_base_range(vcf_line.chr, vcf_line.pos+rep_str.length+1, vcf_line.pos+(rep_str.length*2))
				#STDERR.print "-#{indel_seq}\t#{leading}\t#{rep_str}\t#{trailing}\t"
				if(rep_str == leading || rep_str == trailing)
					vcf_line.info[:"repeat"].push("T")
					vcf_line.info[:"repeatStr"].push(rep_str)
					repeat_count = 0
					repeat_region_length = 0
					#STDERR.print "TRUE\t"
					if(rep_str == leading)
						window = reference.get_base_range(vcf_line.chr, vcf_line.pos - window_size, vcf_line.pos+indel_seq.length)
						if(window =~/(#{rep_str})+$/)
							repeat_region_length = $&.length
							repeat_count = $&.length / rep_str.length
							#STDERR.print "#{$&}\t#{repeat_count}\t"
							#STDERR.print "#{rep_str.length}\t"
							#STDERR.puts "#{repeat_region_length}\t"
						else
							raise "repeat not found when counting repeats!"
						end
					else # trailing
						window = reference.get_base_range(vcf_line.chr, vcf_line.pos+1, vcf_line.pos+1+window_size)
						if(window =~/^(#{rep_str})+/)
							repeat_region_length = $&.length
							repeat_count = $&.length / rep_str.length
							#STDERR.print "#{$&}\t#{repeat_count}\t"
							#STDERR.print "#{rep_str.length}\t"
							#STDERR.puts "#{repeat_region_length}\t"
						else
							raise "repeat not found when counting repeats!"
						end
					end
					vcf_line.info[:"repeatLength"].push(rep_str.length)
					vcf_line.info[:"repeatCnt"].push(repeat_count)
					vcf_line.info[:"repeatRegLength"].push(repeat_region_length)
				else # not repeat
					#STDERR.puts "FALSE"
					vcf_line.info[:"repeat"].push("F")
					vcf_line.info[:"repeatStr"].push('.')
					vcf_line.info[:"repeatCnt"].push('.')
					vcf_line.info[:"repeatLength"].push('.')
					vcf_line.info[:"repeatRegLength"].push('.')
				end
			elsif(alt.length > ref.length)  # insertion
				indel_seq = alt.dup
				indel_seq[0]=''
				rep_str = get_repeat_substr(indel_seq)
				leading = reference.get_base_range(vcf_line.chr, vcf_line.pos-rep_str.length+1, vcf_line.pos)
				trailing = reference.get_base_range(vcf_line.chr, vcf_line.pos+1, vcf_line.pos+rep_str.length)
				#STDERR.print "+#{indel_seq}\t#{leading}\t#{rep_str}\t#{trailing}\t"
				if(rep_str == leading || rep_str == trailing)
					vcf_line.info[:"repeat"].push("T")
					vcf_line.info[:"repeatStr"].push(rep_str)
					repeat_count = 0
					repeat_region_length = 0
					#STDERR.print "TRUE\t"
					if(rep_str == leading)
						window = reference.get_base_range(vcf_line.chr, vcf_line.pos - window_size, vcf_line.pos) + indel_seq
						if(window =~/(#{rep_str})+$/)
							repeat_region_length = $&.length
							repeat_count = $&.length / rep_str.length
							#STDERR.print "#{$&}\t#{repeat_count}\t"
							#STDERR.print "#{rep_str.length}\t"
							#STDERR.puts "#{repeat_region_length}\t"
						else
							raise "repeat not found when counting repeats!"
						end
					else # trailing
						window = indel_seq + reference.get_base_range(vcf_line.chr, vcf_line.pos+1, vcf_line.pos+1+window_size)
						if(window =~/^(#{rep_str})+/)
							repeat_region_length = $&.length
							repeat_count = $&.length / rep_str.length
							#STDERR.print "#{$&}\t#{repeat_count}\t"
							#STDERR.print "#{rep_str.length}\t"
							#STDERR.puts "#{repeat_region_length}\t"
						else
							raise "repeat not found when counting repeats!"
						end
					end
					vcf_line.info[:"repeatLength"].push(rep_str.length)
					vcf_line.info[:"repeatCnt"].push(repeat_count)
					vcf_line.info[:"repeatRegLength"].push(repeat_region_length)
				else
					#STDERR.puts "FALSE"
					vcf_line.info[:"repeat"].push("F")
					vcf_line.info[:"repeatStr"].push('.')
					vcf_line.info[:"repeatCnt"].push('.')
					vcf_line.info[:"repeatLength"].push('.')
					vcf_line.info[:"repeatRegLength"].push('.')
				end
			else
				raise "This variant is not an indel!: #{line}"
			end
		end
		puts vcf_line.print
	end
end
