require 'read'
require 'entropy' 
require 'read_allele'
# this abstract class represents a single indel on a single read from a SAM file
class Read_indel < Read_allele
	attr_reader :length

  
	public

	def initialize( read, offset, start_pos, length )
		super(read, offset, start_pos)
		@length = length
		# @pair_end_ratio = read.i_size()/read.seq.size #TODO support multiple indels
		# in a single read
	end


	def frame_shift
		return 0 if @length%3 == 0
		return 1
	end

	def even_length
		return 1 if @length%2 == 0
		return 0
	end


	def type
		return :indel
	end

	def var_rate
		return (@read.var_count - length).to_f/@read.seq.length.to_f
	end


	# # def in_repeat_region? # 	bed = Bed_file.instance() # 	return 1 if
	# bed.pos_included?(@read.ref, @start_pos) # 	return 0 # end



	# returns the indel's local entropy, checks 5*indel length on both sides,
	# returns the minimum. Uses indel length as entropy window size
	def local_entropy
		seq_length = (5*@length)
		upper_limit = end_pos + seq_length - 1
		ref = Ref_seq.instance
		if(is_deletion)
			lower_limit = @start_pos - seq_length + 1
			lower_limit = 1 if lower_limit < 1
			lower_seq = ref.get_base_range(chr(), lower_limit, end_pos())
			upper_seq = ref.get_base_range(chr(), @start_pos + 1, upper_limit+1)
		else
			lower_limit = @start_pos - seq_length
			lower_limit = 1 if lower_limit < 1
			lower_seq = "#{ref.get_base_range(chr(), lower_limit, @start_pos-1)}#{alt_seq()}"
			upper_seq = "#{alt_seq()[1..-1]}#{ref.get_base_range(chr(), end_pos(), upper_limit)}"
		end
		upstream_entropy = entropy(lower_seq, @length) 
		dnstream_entropy = entropy(upper_seq, @length)
		if upstream_entropy < dnstream_entropy
			return (upstream_entropy*1000).round/1000.0
		else
			return (dnstream_entropy*1000).round/1000.0
		end
	end


	def simple_local_entropy
		ref = Ref_seq.instance
		lower_limit = @start_pos - 10
		lower_limit += 1 if is_deletion
		lower_limit = 1 if lower_limit < 1
		upper_limit = end_pos() + 10
		return (entropy(ref.get_base_range(chr(), lower_limit, upper_limit), @length) * 1000).round/1000.0
	end



	# ## abstract methods ## end_pos read_end_pos indel_bases is_deletion

	private



end #class Read_indel
