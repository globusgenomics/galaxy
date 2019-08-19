class NonRefPosError < Exception
end
class EORefError < Exception
end

# #designed to only load one chromosome into memory at once, chromosome should be given as a symbol
$TEST_REF_PATH = 'test_ref'
class Ref_seq

	$ref_seq_instance = nil

	def initialize(filename, samtools=false)
		@samtools = samtools
		if($DEBUG || samtools)
			@filename = filename
			unless( File.exists?(filename) )
				`samtools faidx #{filename}`
			end
		else
			@bases = Array.new
			@reached_eof = false
			begin
				@file = File.open(filename, 'r')
			rescue
				raise "Could not open reference file #{filename}\n Stopping."
			end

			@line = get_line until @line =~ /^>(\S+)/ || @line == :eof
			# load_next_chr(nil)
		end
		$ref_seq_instance = self
	end


	def self.instance
		raise "ERROR: Trying to access the ref sequence before it is loaded!" if $ref_seq_instance == nil
		return $ref_seq_instance
	end


	# #returns the ascii code of the base
	def get_base(chr, coor)
		if($DEBUG || @samtools)
			# puts "samtools faidx #{@filename} #{chr}:#{coor}-#{coor}"
			return `samtools faidx #{@filename} #{chr}:#{coor}-#{coor}`.split("\n")[1].chomp
		else
			while chr != @chromosome
				load_next_chr(chr)
			end
			if coor < 1 || coor >@bases.length
				STDERR.puts "WARNING: Illegal reference coordinate: #{chr}:#{coor}, returning N" unless $TEST
				return "N"
			end 
			coor = coor - 1 #genome coordinates start counting at 1 instead of 0
			return @bases[coor]
		end
	end

	def get_base_range(chr, start, stop)
		if($DEBUG || @samtools)

			return `samtools faidx #{@filename} #{chr}:#{start}-#{stop}`.split("\n")[1].chomp
		else
			while chr != @chromosome
				load_next_chr(chr)
			end
			if start < 1 || stop > @bases.length
                                STDERR.puts "WARNING: Illegal reference coordinate(s): #{chr}:#{start}-#{stop}"
                                return "N"*(stop-start+1)
                        end
			return @bases[((start-1)..(stop-1))]
		end
	end



	private

	def load_next_chr(chr)
		puts "loading reference chromosome #{chr}...\r"
		if(@line == :eof)
			if(@reached_eof)
				raise EORefError, "Trying to load chromosome #{chr}, but it is not in the reference file!\n"
			end
			@reached_eof = true
			reload()
		end

		@bases = Array.new
                puts "#{@line}"
		if @line =~ /^>(\S+)/
			@chromosome = $1.sub("chr","").to_sym
			@bases = ""
		else
			raise "Missing header line in reference file, instead found: #{@line}"
		end
		@line = get_line
		if( chr != nil && @chromosome != chr )
			@line = get_line until @line =~ /^>(\S+)/ || @line == :eof
			return load_next_chr(chr)
		end
		# found the right chromosome
		@reached_eof = false
		STDOUT.flush
		while( @line != :eof )
			if @line =~ /^>(\S+)/
				break
			else
				@bases << @line
			end
			@line = get_line
		end
		@bases.upcase!
		STDOUT.flush
	end



	def reload
		@file.rewind
		@line = get_line until @line =~ /^>(\S+)/ || @line == :eof
	end


	def get_line
		begin
			line = @file.readline
			line.chomp!
		rescue EOFError
			line = :eof
		end
		return line
	end

end
