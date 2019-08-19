class NonRefPosError < Exception
end
class EORefError < Exception
end
class FaiItem < Struct.new(:chrom, :length, :offset, :bases_per_line, :bytes_per_line)
end

# #designed to only load one chromosome into memory at once, chromosome should be given as a symbol
$TEST_REF_PATH = 'test_ref'
class Ref_seq
	attr_reader :bases
	$ref_seq_instance = nil

	def initialize(filename)
		#@samtools = samtools
		if($DEBUG )
			@filename = filename
			unless( File.exists?(filename) )
				`samtools faidx #{filename}`
			end
		else
			@fai = nil
			if( File.exists?(filename + '.fai') )
				@fai = load_fasta_index(filename + '.fai')
			end
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
		raise "ERROR: Trying to access the ref sequence before it is loaded!" if $ref_seq_instance.nil?
		return $ref_seq_instance
	end


	# #returns the ascii code of the base
	def get_base(chr, coor)
		if($DEBUG)
			# puts "samtools faidx #{@filename} #{chr}:#{coor}-#{coor}"
			x = `samtools faidx #{@filename} #{chr}:#{coor}-#{coor}`.split("\n")[1]
			if(x.nil?)
				# raise NonRefPosError, "Illegal reference coordinate: #{chr}:#{coor}" 
				return "N"
				
			else
				return x.chomp
			end
		elsif( @fai != nil )
			index = @fai[chr]
			coor -= 1
			@file.seek(index.offset + coor / index.bases_per_line * index.bytes_per_line + coor % index.bases_per_line)
			base = @file.read(1)
			if(base.nil?)
				STDERR.puts "WARNING: Illegal reference coordinate: #{chr}:#{coor}, returning N" unless $TEST
				return "N"
			end
			return base
		else
			while chr != @chromosome
				load_next_chr(chr)
			end
			if coor < 1 || coor >@bases.length
				STDERR.puts "WARNING: Illegal reference coordinate: #{chr}:#{coor}, returning N" unless $TEST
				return "N"
			end 
			return @bases[coor-1]
		end
	end

	def get_base_range(chr, start, stop)
		if($DEBUG )

			return `samtools faidx #{@filename} #{chr}:#{start}-#{stop}`.split("\n")[1].chomp
		elsif( @fai != nil )
			index = @fai[chr]
			start -= 1 # convert to 0-based coordinates
			stop -= 1
			@file.seek(index.offset + start / index.bases_per_line * index.bytes_per_line + start % index.bases_per_line)
			output = ""
			while output.length <= stop-start
				nextline = get_line
				if( nextline =~ /^>(\S+)/ || nextline == :eof )
					STDERR.puts "WARNING: Illegal reference coordinate: #{chr}:#{start+output.length}, returning N" unless $TEST
					while output.length <= stop-start
						output << "N"
					end
				else
					output << nextline
				end
			end
			return output[0..(stop-start)]
		else
			while chr != @chromosome
				load_next_chr(chr)
			end
			if stop < 1 || start > @bases.length
                                STDERR.puts "WARNING: Illegal reference coordinate(s): #{chr}:#{start}-#{stop}"
                                return "N"*(stop-start+1)
			elsif start < 1 && stop > 0
                                STDERR.puts "WARNING: Illegal reference coordinate(s): #{chr}:#{start}-#{stop}"
				return "N"*(start.abs + 1) + @bases[(1..(stop-1))]
			elsif stop > @bases.length 
                                STDERR.puts "WARNING: Illegal reference coordinate(s): #{chr}:#{start}-#{stop}"
				return @bases[((start-1)..(@bases.length-1))] + "N"*(stop - @bases.length)
                        end

			return @bases[((start-1)..(stop-1))]
		end
	end



	private

	def load_fasta_index(faifilename)
		index = Hash.new
		File.foreach(faifilename) do |line|
			match = /(.+?)\t(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/.match(line)
			if( match != nil ) 
				chrom, len, offset, bases_len, bytes_len = match[1,5]
				chrom = chrom.sub("chr","")
				index[chrom] = FaiItem.new(chrom, len.to_i, offset.to_i, bases_len.to_i, bytes_len.to_i)
			else
				STDERR.puts "WARNING: Ignoring unrecognized line in fasta index (.fai) file:\n#{line}"
			end
		end
		return index
	end

	def load_next_chr(chr)
		print "loading reference chromosome #{chr}...                                                                                    \r"
		if(@line == :eof)
			if(@reached_eof)
				raise EORefError, "Trying to load chromosome #{chr}, but it is not in the reference file!\n"
			end
			@reached_eof = true
			reload()
		end

		@bases = Array.new
		if @line =~ /^>(\S+)/
			@chromosome = $1.sub("chr","")
			@bases = ""
		else
			raise "Missing header line in reference file, instead found: #{@line}\nchr=#{chr}, @chromosome=#{@chromosome}"
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
		@line = get_line
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
