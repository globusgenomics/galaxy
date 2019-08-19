class NoCigarError < Exception
end
class FormatError < Exception
end
class FlaggedError < Exception
end

require 'ref_seq'
require 'entropy'

# A single sequencing read, designed to work with the SAM format
class Read
	@@flag_filt = 1796


	attr_accessor :name, :flag, :ref, :pos, :map_qual, :cigar, :i_size,
	  :seq, :qual, :score, :dir
	attr_reader :is_clipped, :color_correction, :print_chr, :var_count

	def initialize( sam_line, use_orig_base_qual, remove_clips=true)
		sam_line.chomp!
		@remove_clips = remove_clips
		read_sam_cols = sam_line.split
		unless read_sam_cols.length > 10
			raise FormatError.new("Malformed SAM line:\nskipping...")
		end
		@name = read_sam_cols[0]
		@flag = read_sam_cols[1].to_i
		if(flag & @@flag_filt != 0)
			raise FlaggedError.new( "Read fails bitwise flag mask" )
		end
		if(read_sam_cols[2].to_s[0,3] == 'chr') #preappend the 'chr' label to the chromosomes in the final output
			@print_chr = true
		else
			@print_chr = false
		end
		@ref = read_sam_cols[2].sub("chr","").to_sym
		@pos = read_sam_cols[3].to_i
		@map_qual = read_sam_cols[4].to_i
		if(sam_line =~ /CM:i:(\d+)/)
			@color_correction = $1.to_i
		else
			@color_correction = 0
		end
		
		# @mate_ref = ref_sam_cols[6] @mate_pos = ref_sam_cols[7] @i_size =
		# read_sam_cols[8]
		@seq = read_sam_cols[9]
		qual = read_sam_cols[10]
		@score = 'NA'
		# parse optional tags of SAM
		orig_qual = nil
		num_fields = read_sam_cols.length
		(11..(num_fields-1)).each do |col_num|
			if read_sam_cols[col_num] =~ /^AS:i:(\S+)$/
				@score = $1
			end
			if read_sam_cols[col_num] =~ /^OQ:Z:([^\t]*)/
				orig_qual = $1
			end
			if(read_sam_cols[col_num] =~ /^NM:i:(\S+)$/)
				@var_count = $1.to_i
			end
		end
		begin
			if(use_orig_base_qual && orig_qual != nil)
				parse_qual(orig_qual)
			else
				parse_qual(qual)
			end
		rescue
			raise "Error parsing read quality. Line:\n#{sam_line}\n\n#{$!}"
		end
		begin
			parse_cigar(read_sam_cols[5])
		rescue NoCigarError
			STDERR.puts "Read #{@name}, #{@ref}:#{@pos} has no cigar code, skipping...\n"
			@cigar = Array.new
		end
		@var_count = var_count() if @var_count.nil?
		# parse bitwise flag to get strand dir
		if @flag.to_s(2).reverse[4,1] == '1'
			@dir = :'-'
		else
			@dir = :'+'
		end

		raise "Illegal SAM line, quality entry length #{@qual_array.length} does not match sequence length #{@seq.length}" if @qual_array.length != @seq.length

	end  #initialize

  
	# returns an array of numerical quality values parsed from the ascii values
	def qual_array
		return @qual_array
	end

	#this may only be run AFTER get_indels is run
	def end_pos
		return @end_pos
	end

	def mean_base_qual
		return (@qual_array.inject(:+).to_f / @qual_array.length.to_f).to_i
	end

	def length
		return @seq.length
	end



	# returns the number of variations/the read length, also adds counts to the variant_counter
	def var_rate()
		nearby_var_sites = Array.new
		count = 0
		read_pos = 0
		ref_pos = 0
		@cigar.each do |entry|
			length = entry[0].to_i
			op = entry[1]
			if(op == :"I")
				nearby_var_sites.push("#{@ref}:#{@pos+ref_pos-1}I#{seq[read_pos..(read_pos+length-1)]}")
				count += 1
				read_pos += length
			elsif(op == :"D")
				nearby_var_sites.push("#{@ref}:#{@pos+ref_pos-1}D#{length}")
				count += 1
				ref_pos += length
			elsif(op == :"P")
				ref_pos += length
			elsif(op == :"M" || op == :"=" || op == :"X")
				(read_pos...(read_pos+length)).each do |j| #count SNPs
					#puts "read_base:#{@seq[j].chr}  ref_base:#{Ref_seq.instance.get_base(@ref, @pos+ref_pos)}"
					begin
						if(op == :"X" ||  (op != :"=" && @seq[j] != Ref_seq.instance.get_base(@ref, @pos+ref_pos)))
							nearby_var_sites.push("#{ref}:#{@pos+ref_pos}#{@seq[j]}")
							count += 1
						end
					rescue NonRefPosError
						# read is outside of the reference genome
						count += 1
					end
					ref_pos += 1
				end
				read_pos += length
			elsif(op == :"N")
				(read_pos...(read_pos+length)).each do |j| 
					ref_pos += 1
				end
				read_pos += length
			else
				raise "Unsupported cigar entry #{op}"
			end
		end
		return ((count.to_f/@seq.length.to_f)*10000).round/10000.0, nearby_var_sites
	end



	def var_count
		count = 0
		read_pos = 0
		ref_pos = 0
		@cigar.each do |entry|
			length = entry[0].to_i
			op = entry[1]
			if(op == :"I")
				count += length
				read_pos += length
			elsif(op == :"D")
				count += length
				ref_pos += length
			elsif(op == :"P")
				ref_pos += length
			elsif(op == :"M" || op == :"=" || op == :"X")
				(read_pos...(read_pos+length)).each do |j| #count SNPs
					#puts "read_base:#{@seq[j].chr}  ref_base:#{Ref_seq.instance.get_base(@ref, @pos+ref_pos)}"
					begin
						if(op == :"X" ||  (op != :"=" && @seq[j] != Ref_seq.instance.get_base(@ref, @pos+ref_pos)))
							count += 1
						end
					rescue NonRefPosError
						# read is outside of the reference genome
						count += 1
					end
					ref_pos += 1
				end
				read_pos += length
			elsif(op == :"N")
				(read_pos...(read_pos+length)).each do |j| #count SNPs
					ref_pos += 1
				end
				read_pos += length
			else
				raise "Unsupported cigar entry #{op}"
			end
		end
		return count.to_f
	end


	# returns the number of gaps/the read length
	def gap_rate
		count = 0
		@cigar.each do |entry|
			op = entry[1]
			if(op == :"I")
				count += 1
			elsif(op == :"D")
				count += 1
			elsif(op == :"M" || op == :"N" || op == :"=" || op == :"X" || op == :"P")
				#do nothing
			else
				raise "Unsupported cigar entry #{op}"
			end
		end
		return ((count.to_f/@seq.length.to_f)*10000).round/10000.0
	end



	def get_indels
		indels = Array.new
		offset = 0
		pos = @pos
		@cigar.each do |entry|
			begin
				if(entry[1]==:"D")
					indels.push(Read_deletion.new(self, offset-1, pos-1, entry[0]))
					pos += entry[0]
				elsif(entry[1]==:"P")
					pos += entry[0]
				elsif(entry[1]==:"I")
					indels.push(Read_insertion.new(self, offset-1, pos-1, entry[0]))
					offset += entry[0]
				elsif(entry[1]==:"M" || entry[1] == :"N" || entry[1] == :"=" || entry[1] == :"X")
					offset += entry[0]
					pos += entry[0]
				else
					raise "Unsupported cigar entry: #{entry}"
				end
			rescue NonRefPosError
				# indel is outside of the reference genome, and therefore not reliable so skip it
			end
		end
		@end_pos = pos - 1
		return indels
	end
	

	# returns the indels in the read, also updates the depth coverage Hashes. NOTE: depths are for indels, are NOT applicable to SNPs # this is probably the slowest method, needs to be optimized, especially the range.each sections
	def get_indels_with_depth(depth, ref_depth)
		indels = Array.new
		offset = 0
		pos = @pos
		prev_ins = false # was the previous cigar entry an insertion?
		first_read_base = true
		@cigar.each do |entry|
			begin
				if(entry[1]==:"M" || entry[1]==:"N" || entry[1]==:"=" || entry[1]==:"X")
					((pos-1)..(pos+entry[0]-2)).each do |i|
						if(first_read_base)# skip the first base of the read. This is required for consistant insertion depths
							first_read_base = false
							next
						end
						if(depth[i].nil?)
							ref_depth[i] = 1  
							depth[i] = 1
						else
							ref_depth[i] += 1
							depth[i] += 1
						end
						if(prev_ins) # if the previus entry was an insertion, this site has already been counted
							ref_depth[i] -= 1 
							depth[i] -= 1
							prev_ins = false
						end
					end
					offset += entry[0]
					pos += entry[0]
				else
					if(entry[1]==:"D")
						((pos-1)..(pos+entry[0]-2)).each do |ipos| # count total depth across deletions
							if(depth[ipos].nil?)
								depth[ipos]=1
								ref_depth[ipos]=0
							else
								depth[ipos] += 1
							end
						end
						indels.push(Read_deletion.new(self, offset-1, pos-1, entry[0]))
						pos += entry[0]
						prev_ins = false
					elsif(entry[1]==:"P")
						pos += entry[0]
						prev_ins = false
					elsif(entry[1]==:"I")
						if(depth[pos-1].nil?)
							ref_depth[pos-1] = 0
							depth[pos-1]=1
						else
							depth[pos-1]+=1
						end
						indels.push(Read_insertion.new(self, offset-1, pos-1, entry[0]))
						offset += entry[0]
						prev_ins = true
						first_read_base=false
					else
						raise "Unsupported cigar entry: #{entry}"
					end
				end
			rescue NonRefPosError
				# indel is outside of the reference genome, and therefore not reliable so skip it
				prev_ins = false
			end
		end
		@end_pos = pos - 1
		return indels
	end
	

# return the fraction of the read sequence that is G or C
	def gc_content
		return ((@seq.count('GCgc').to_f/@seq.length.to_f)*10000).round/10000.0
	end


	def base_qual_sum
		return @qual_array.inject(:+)
	end


	def complexity
		#return 0.5
		chunk = nil
		max_rep_score = 0
		# max_chunk = 0
		max = 16
		max = @seq.length/2 if max > @seq.length/2
		(2..max).each do |window|
			(0..(@seq.length-window)).each do |i|
				#rep_score = 0
				chunk = @seq[i,window]
				#rep_score = @seq.scan(chunk).length * window
				index = -window
				count = 0
				while (index = @seq.index(chunk, index+window))
					count += 1
				end
				rep_score = count * window
				# # 	max_chunk = chunk if rep_score > max_rep_score && rep_score > window
				max_rep_score = rep_score if rep_score > max_rep_score && rep_score > window
			end
		end
		max_rep_score = max_rep_score.to_f/@seq.length.to_f
		raise "repetiveness should not be over 1.  #{max_rep_score}  -- #{@seq}" if max_rep_score > 1
		# return "#{max_chunk}: #{max_rep_score}"
		return max_rep_score
	end

	def read_entropy(window_size)
		return (entropy(@seq, window_size)*1000).round/1000.0
	end

	def normalized_entropy(window_size)
		return (entropy(@seq, window_size)/(window_size * Math.log(4)) *1000).round/1000.0
	end


	private #################
	

	def parse_qual(qual)
		@qual_array = []
		qual.each_byte do |ascii|
			@qual_array.push(ascii - 33)
		end
	end


	def parse_cigar(cigar_str)
		@is_clipped = false
		raise NoCigarError, "no cigar code" if cigar_str == nil || cigar_str == "" || cigar_str == "*"
		@cigar=Array.new
		sum = 0
		cigar_str.scan(/(\d+)([MIDNSHPX=])/) do |entry|
			if( @remove_clips )
				if(entry[1] == "S") # soft clip
					@is_clipped = true
					if(@cigar.length == 0) # clip is at beginning
						@seq.slice!(0,entry[0].to_i) #cut the clipped sequence out of the read
						@qual_array.slice!(0,entry[0].to_i)
					else # clip is at the read's end
						@seq.slice!(-(entry[0].to_i),entry[0].to_i) #cut the clipped sequence out of the read if entry[1] == "M" || entry[1] == "I"
						@qual_array.slice!(-entry[0].to_i,entry[0].to_i)
					end
					next
				elsif(entry[1] == 'H') # hard clip
					@is_clipped = true
					next # just skip hard clips
				end
			end
			sum += entry[0].to_i if entry[1] == "M" || entry[1] == "I" || entry[1] == "N" || entry[1]=="X" || entry[1]=="="
			@cigar.push([entry[0].to_i, entry[1].to_sym])
		end
		raise "Cigar does not match sequence length! #{sum} != #{@seq.length}" if sum != @seq.length
	end



	
  


end
