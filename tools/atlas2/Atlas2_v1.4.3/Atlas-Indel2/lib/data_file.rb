#!/usr/bin/ruby
#methods to open and process general data files
#usage: require 'Data_file'  --   file=Data_file.new(file.vcf, output.vcf, ",")  --   file.process_file { |cols, writer, line| writer.puts cols.join(",") }

class Data_file

	def initialize(filename, outfile=nil, sep="\t")
		@file = File.open(filename, 'r')
		if(outfile==nil)
			@writer=STDOUT
		else
			@writer=File.open(outfile, 'w')
		end
		@sep=sep
	end

	#method that processes the file.  Take a code block that is responsible for writing output
	def process(&block)
		@file.each_line do |line|
			line.chomp!
			cols = line.split(@sep)
			if block_given?
				case block.arity
					when 1
						yield cols
					when 2
						yield cols, @writer
					when 3
						yield cols, @writer, line
				end	
			end
		end
	end

	def close
		@file.close
		unless @writer.tty? #can't close a terminal
			@writer.close
		end
	end
end
