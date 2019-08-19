require 'read'
require 'entropy' 
# this abstract class represents a single allele (ref/SNP/indel) on a single read from a SAM file
class Read_allele
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


	# returns the distance to the next polymorphism
	def proximity
		cigar = @read.cigar
		cigar_i = nil
		read_pos = 0
		(0...(cigar.length)).each do |i|
			if(read_pos == @offset + 1)
				cigar_i = i
				break
			end
			read_pos += cigar[i][0] unless cigar[i][1] == :"D"
		end

		lower_dist, found_lower_morph = get_morph_dist(cigar, true, cigar_i)
		upper_dist, found_upper_morph = get_morph_dist(cigar, false, cigar_i)

		if(found_lower_morph)
			if(lower_dist < upper_dist || !found_upper_morph)
				return lower_dist
			else
				return upper_dist
			end
		elsif(found_upper_morph)
			return upper_dist
		end

		# did not find any other variations on this read
		return lower_dist if lower_dist > upper_dist
		return upper_dist
	end


	def dist_3
		return @offset + 1 if @read.dir == :'-'
		return @read.length - @offset -1
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

	def nqs #require quality >= 15 for surrounding 5 bases
		#check 5 upstream bases
		lower = @offset-4
		upper = @offset
		if(upper >= 0)
			lower = 0 if lower < 0
			begin
				@read.qual_array[lower..upper].each do |qual|
					return 0 if qual < 15
				end
			rescue
				raise "Error with nqs check: lower:#{lower} upper:#{upper} qual_array:#{@read.qual_array}  --  #{@read.qual_array[lower..upper]}\n#{$!}"
			end
		end
		#check 5 downstream bases
		lower = @offset+1
		upper = @offset +5
		if(lower < @read.qual_array.length)
			upper = @read.qual_array.length - 1  if upper >= @read.qual_array.length
			@read.qual_array[lower..upper].each do |qual|
				return 0 if qual < 15
			end
		end
		return 1
	end


    def avnbq
        #Returns average base quality score around the allele
        lower = @offset-4
        upper = @offset+4
        if lower < 0
            lower=0
        end
        if upper > @read.qual_array.length - 1
            upper=@read.qual_array.length - 1
        end
        total_qual=0
        @read.qual_array[lower..upper].each do |qual|
            total_qual = total_qual+qual
        end
        return total_qual/(upper-lower+1)
    end

	private


	def get_morph_dist(cigar, is_lower, cigar_i)
		found = false
		dist = 0
		if(is_lower)
			dir = -1
		else
			dir = 1
		end
		entry = cigar[cigar_i+(1*dir)]
		if(entry != nil && entry[1] == :"M")
			ref = Ref_seq.instance
			(0...(entry[0])).each do |i|
				if(is_lower)
					start = @start_pos
					start += 1 if is_deletion
					read_start = read_start_pos
				else
					start = end_pos
					read_start = read_end_pos
				end
				start += (1*dir) if cigar[cigar_i][1] == :"D" #adjust ref_index by 1 for deletions
				# puts "ref[#{start+(i*dir)}] #{ref.get_base(@read.ref, start+(i*dir))}
				# read[#{read_start+(i*dir)}] #{@read.seq[read_start+(i*dir), 1]}"
				if(ref.get_base(@read.ref, start+(i*dir)) != @read.seq[read_start+(i*dir)])
					found = true
					break
				end
				dist += 1
			end
		end
		if(!found && cigar_i+(2*dir) > 0 && cigar_i+(2*dir) < cigar.length)
			found = true
		end
		return dist, found
	end

end #class read_variant
