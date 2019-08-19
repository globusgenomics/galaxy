require 'read'
require 'entropy' 
# this abstract class represents a single variant on a single read from a SAM file
class Read_variant
	attr_reader :length, :start_pos, :read

  
	public

	def initialize( read, offset, start_pos )
		@read = read
		@offset = offset # how far from the start of the read
		@start_pos = start_pos  # the coordinate of the variant in the genome
		# @pair_end_ratio = read.i_size()/read.seq.size #TODO support multiple indels
		# in a single read
	end


	def map_qual
		return @read.map_qual
	end

	def chr
		return @read.ref
	end


	def read_start_pos
		return @offset
	end


	# # def in_repeat_region? # 	bed = Bed_file.instance() # 	return 1 if
	# bed.pos_included?(@read.ref, @start_pos) # 	return 0 # end


	# returns 1 if the variant is within 5bp of a read end
	def near_read_end
		if(@offset < 4 || @read.seq.length - read_end_pos < 5)
			return true
		end
		return false
	end


	private



end #class read_variant
