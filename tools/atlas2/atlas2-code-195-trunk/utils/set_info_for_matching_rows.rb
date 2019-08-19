#!/usr/bin/ruby
#add/change an info entry for lines that exactly match in Locus and allele.  Doesn't work for multi-allelic files

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Set_info

	def initialize(file, file2, info_tag, info_value)
		@filename = file
		@filename2 = file2
		@source_header = `grep "#CHROM	" #{file2}`.chomp.split("\t")
		@info_tag = info_tag.to_sym
		@info_value = info_value
	end

	def run
		labels=nil
		File.open(@filename, 'r').each_line do |line|
			line.chomp!
			if(line[0,1] == '#')
				if( line[0,4] == '#CHR' )
					labels = line.split("\t") 
				end
				puts line
				next
			end
			vcf_line = Vcf_line.read_line(line, false, labels)
			edit_line(vcf_line)
			puts vcf_line.print()
		end
	end

	def edit_line(line)
		source_str=`grep -P "^#{line.chr}	#{line.pos}	.	#{line.ref}	#{line.alt[0]}	" #{@filename2}`.chomp
		if(source_str.nil? || source_str == "")
			STDERR.puts "WARNING: could not find line #{line.chr}:#{line.pos} in the source VCF"
			return
		end
		line.info[@info_tag]=@info_value
	end
end


if(ARGV.length < 3)
	puts "USAGE: ruby set_info_for_matching_rows.rb vcf_to_edit source_vcf info_tag info_value"
	exit 1
end

editor = Set_info.new(ARGV[0], ARGV[1], ARGV[2], ARGV[3])
editor.run
