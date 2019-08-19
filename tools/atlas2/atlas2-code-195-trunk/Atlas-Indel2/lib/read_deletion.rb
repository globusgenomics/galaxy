require 'entropy'

class Read_deletion < Read_indel
	def initialize( read, offset, start, length )
		super(read, offset, start, length)
	end




	def neighbor_entropy
		bases = Ref_seq.instance.get_base_range(@read.ref, @start_pos-10, end_pos+10)
		return entropy(bases, @length)
	end



	def ref_seq
		return Ref_seq.instance.get_base_range(@read.ref, @start_pos, end_pos)
	end


	def alt_seq
		return Ref_seq.instance.get_base(@read.ref, @start_pos)
	end

	def end_pos
		return @start_pos+@length
	end


	def read_end_pos
		return @offset + 1
	end


	def is_deletion
		return true
	end


	def possible_homopolymer
		del_bases = ref_seq
		str = "#{@read.seq[@offset,1]}#{del_bases}#{@read.seq[read_end_pos,1]}"
		str.squeeze!
		return 0 if str.length > 2
		return 1
	end



	def to_s
		return "#{@read.ref}:#{@start_pos}D#{@length}"
	end

    
end
