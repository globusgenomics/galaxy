require 'vcf_line.rb'
class Vcf_file
	attr_reader :labels, :sample_names, :filename

	def initialize(filename=nil)
		@filename = filename
		if(filename.nil?)
			@file = STDIN
		else
			@file = File.open(filename, 'r')
		end
		@header = Array.new
		load_header()
	end

	def each_line
		while(!eof())
			yield getline()
		end
	end

	def getline()
		vcf_str = @file.gets
		return Vcf_line.read_line(vcf_str, false, @labels)
	end


	def print_header
		return @header.join("\n") + "\n" + @labels.join("\t")
	end

	def eof()
		return @file.eof()
	end

	def close()
		return @file.close()
	end

	private

	def load_header
		begin
			line = @file.gets
			while(line[0,1] == '#')
				line.chomp!
				if( line[0,4] == '#CHR' )
					@labels = line.split("\t") 
					return
				else
					@header.push(line)
					line = @file.gets
				end
			end
		rescue
			raise "No valid header in the VCF file #{@filename}"
		end
	end
end
