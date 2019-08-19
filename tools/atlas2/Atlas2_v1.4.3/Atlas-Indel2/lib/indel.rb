# this class represents a single indel as determined from one of many read_indels

class Indel
	attr_reader :read_count, :num_passed_reads, :sample, :max_map_qual, :nearby_var_counts
	attr_accessor :req_include
	def initialize(read_indel, z_score, cutoff, sample=nil)
		if(read_indel.read.dir == :'-')
			@plus_strand_dir = false
			@minus_strand_dir = true
		else
			@plus_strand_dir = true
			@minus_strand_dir = false
		end
		if(read_indel.near_read_end)
			@near_read_end_count=1
		else
			@near_read_end_count=0
		end
		#		@z_score = z_score
		@indel = read_indel
		@read_count = 1
		@map_qual = read_indel.map_qual
		@ave_nqs = read_indel.avnbq
		@var_rate, nearby_var_sites = read_indel.read.var_rate
		@var_rate2 = read_indel.var_rate
		@cutoff = cutoff
		@sample = sample
		@color_correction = read_indel.read.color_correction
		@max_map_qual = read_indel.map_qual
		@high_var_count = 0
		@high_var_count += 1 if read_indel.read.var_count - read_indel.length > 1
		@req_include = false
		@nearby_var_counts = Hash.new
		# raise "failed assertion!" if !nearby_var_sites.include?(read_indel.to_s)
		nearby_var_sites.delete(read_indel.to_s) # remove the variant site
		nearby_var_sites.each {|pos| @nearby_var_counts[pos]=1 unless pos == @indel.start_pos}
		#		if(z_score > cutoff)
		#			@num_passed_reads = 1
		#		else
		#			@num_passed_reads = 0
		#		end
	end

	# precond: the added read_indel == the current read_indel (but are different reads)
	def  add_read_indel(z_score, read_indel)
		if(read_indel.read.dir == :'-')
			@minus_strand_dir = true
		else
			@plus_strand_dir = true
		end
		if(read_indel.near_read_end)
			@near_read_end_count+=1
		end
		@max_map_qual = read_indel.map_qual if read_indel.map_qual > @max_map_qual
		#		@z_score += z_score
		@read_count += 1
		@map_qual += read_indel.map_qual
		@ave_nqs += read_indel.avnbq
		v, nearby_var_sites = read_indel.read.var_rate
		@var_rate += v
		@var_rate2 += read_indel.var_rate
		@color_correction += read_indel.read.color_correction
		@high_var_count += 1 if read_indel.read.var_count - read_indel.length > 1
		# raise "failed assertion!" if !nearby_var_sites.include?(read_indel.to_s)
		nearby_var_sites.delete(read_indel.to_s) # remove the variant sites
		nearby_var_sites.each do |pos|
			if(@nearby_var_counts[pos].nil?)
				@nearby_var_counts[pos]=1 unless pos == @indel.start_pos
			else
				@nearby_var_counts[pos]+=1 unless pos == @indel.start_pos
			end
		end
		#		if(z_score > @cutoff)
		#			@num_passed_reads += 1
		#		end
	end

	def high_var_ratio
		return ((@high_var_count.to_f / @read_count.to_f)*10000.0).round/10000.0
	end
	
	def artifact_variant_ratio
		return @nearby_var_counts.values.max.to_f / @read_count.to_f
	end
	
	def artifact_variants
		i = 0
		@nearby_var_counts.each_value do |count|
			i+=1 if count.to_f / @read_count.to_f > 0.3
		end
		return i
	end

	def z_score
		return ((@z_score.to_f / @read_count.to_f)*10000.0).round/10000.0
	end

	def sum_z_scores
		return (((@z_score.to_f*100.0).round/100.0)*10000.0).round/10000.0
	end

	def read_indel
		return @indel
	end

	def all_near_read_end
		return 1.0 if @near_read_end_count == @read_count
		return 0.0
	end

	def near_read_end_ratio
		return @near_read_end_count.to_f / @read_count.to_f
	end

	def near_read_end_and_strand_bais
		if(all_near_read_end == 1.0 && fails_strand_dir_filt)
			return 1.0
		end
		return 0.0
	end

	def mean_color_corrections
		return @color_correction.to_f/@read_count.to_f
	end

	def mean_map_qual
		return ((@map_qual.to_f / @read_count.to_f)*100.0).round/100.0
	end

	def mean_avnqs
		return ((@ave_nqs.to_f / @read_count.to_f)*100.0).round/100.0
	end

	def length
		return @indel.length
	end

	def frameshift
		return @indel.frame_shift
	end

	def mean_var_rate
		return ((@var_rate / @read_count.to_f)*100.0).round/100.0
	end

	def mean_var_rate2
		return ((@var_rate2 / @read_count.to_f)*100.0).round/100.0
	end

	def simple_local_entropy
		return @indel.simple_local_entropy
	end

	def local_entropy
		return @indel.local_entropy
	end

	# return true if variant reads are only in one strand direction
	def fails_strand_dir_filt
		return !(@plus_strand_dir && @minus_strand_dir)
	end

	def new_fails_strand_dir_filt
		return false if @read_count < 5
		return !(@plus_strand_dir && @minus_strand_dir)
	end

	# returns the 0 if there is a read in each strand direction, and the variant read depth otherwise (higher is worse)
	def strand_score
		if((@plus_strand_dir && @minus_strand_dir))
			return 0.0
		else
			return @read_count.to_f
		end
	end

	def strand_dir
		return 0.0 if fails_strand_dir_filt
		return 1.0
	end

	def new_strand_dir
		return 0.0 if new_fails_strand_dir_filt
		return 1.0
	end

	def site_passed_over_var
		return @num_passed_reads.to_f / @read_count.to_f
	end

	def is_deletion
		return @indel.is_deletion
	end

	# other must be a read_indel
	#	def == (other)
	#		if(@indel.is_deletion && other.is_deletion)
	#			if(@indel.start_pos == other.start_pos && @indel.ref_seq == other.ref_seq)
	#				return true
	#			end
	#		elsif( ! @indel.is_deletion && ! other.is_deletion)
	#			if(@indel.start_pos == other.start_pos && @indel.alt_seq == other.alt_seq)
	#				return true
	#			end
	#		else
	#				return false
	#		end
	#	end

end

