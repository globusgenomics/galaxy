require 'read_indel'
require 'assertion_failure'

class Read_insertion < Read_indel
	attr_reader :alt_seq

	def initialize( read, offset, start, length )
		super(read, offset, start, length)
		ref = Ref_seq.instance
		begin
			@alt_seq = "#{ref.get_base(read.ref, start)}#{@read.seq[@offset+1,@length]}"
		end
	end


	def ref_seq
		return Ref_seq.instance.get_base(@read.ref, @start_pos)
	end



	def dist_3
		return @offset + 1 if @read.dir == :'-'
		return @read.length - (@offset + 1 +@length)
	end



	def base_qual
		begin
			return @read.qual_array[@offset+1] if @length == 1
			quals = @read.qual_array[(@offset+1)..(@offset+@length)]
			return (quals.inject(0.0) { |sum, el| sum + el } / quals.size).round
		rescue
			raise "Error calculating base_qual #{@offset}...#{@offset+@length-1}  length:#{@length} #{quals}"
		end
	end




	def nqs #require quality >= 15 for surrounding 5 bases
		# #check 5 upstream bases
		lower = @offset-4
		upper = @offset
		if(upper >= 0)
			lower = 0 if lower < 0
			raise "problem calculating nqs, lower:#{lower} upper:#{upper}" if (lower..upper) == nil
			@read.qual_array[(lower..upper)].each do |qual|
				return 0 if qual < 15
			end
		end
		# #check 5 downstream bases
		lower = @offset+@length+1
		upper = @offset +@length+5
		if(lower < @read.qual_array.length)
			upper = @read.qual_array.length - 1  if upper >= @read.qual_array.length
			raise "problem calculating nqs, lower:#{lower} upper:#{upper}" if (lower...upper) == nil
			@read.qual_array[lower..upper].each do |qual|
				return 0 if qual < 15
			end
		end
		return 1
	end


	def avnbq
		# Returns average base quality score around the indel
		lower = @offset-4
		upper = @offset+4+@length
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

    
	def end_pos
		return @start_pos+1
	end


	def is_deletion
		return false
	end


	def read_end_pos
		return @offset + @length + 1
	end


	def possible_homopolymer
		if(@offset < 0)
			start=Ref_seq.instance.get_base(@read.ref, @start_pos)
			str = "#{start}#{@read.seq[0..read_end_pos]}"
		elsif(@offset + @length + 1 >= @read.length)
			tail=Ref_seq.instance.get_base(@read.ref, end_pos)
			str = "#{@read.seq[(@offset)..read_end_pos]}#{tail}"
		else
			str = @read.seq[(@offset)..read_end_pos]
		end
		assert(str.length == 2+@length)
		str.squeeze!
		return 0 if str.length > 2
		return 1
	end


	def to_s
		return "#{@read.ref}:#{@start_pos}I#{@alt_seq[1..(@alt_seq.length)]}"
	end

end
