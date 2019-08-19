#!/usr/bin/ruby

$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
class Vcf_editor

	def initialize(file=nil)
		@filename = file
	end

	def run
		labels=nil
		if(@filename.nil?)
			reader = STDIN
		else
			reader= File.open(@filename, 'r')
		end
		reader.each_line do |line|
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
		#abstract method, do nothing
	end
end
if __FILE__ == $0
	if(ARGV.length < 1)
		puts "USAGE: ruby copy_editor.rb vcf_file"
		exit 1
	end
	
	editor = Vcf_editor.new(ARGV[0])
	editor.run
end
