#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Copy_info

	def initialize(file2, info_tag)
		@filename2 = file2
		@source_header = `grep "#CHROM	" #{file2}`.chomp.split("\t")
		@info_tags = info_tag.split(",")
		@info_tags.map! {|x| x.to_sym}
	end

	def run
		labels=nil
		STDIN.each_line do |line|
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
		source_str=`grep -P "^#{line.chr}	#{line.pos}	[^	]*	#{line.ref}	#{line.alt[0]}	" #{@filename2}`.chomp
		if(source_str.nil? || source_str == "")
			STDERR.puts "WARNING: could not find line #{line.chr}:#{line.pos} in the source VCF"
			return
		end
		source_line = Vcf_line.read_line(source_str, false, @source_header)
		@info_tags.each do |tag|
			line.info[tag]=source_line.info[tag]
		end
	end
end


if(ARGV.length < 2 || ARGV[0] == "-h" || ARGV[0] == "--help")
	puts "USAGE: cat vcf_to_edit | ruby copy_info.rb source_vcf info_tags(comma separated)"
	exit 1
end

editor = Copy_info.new(ARGV[0], ARGV[1])
editor.run
