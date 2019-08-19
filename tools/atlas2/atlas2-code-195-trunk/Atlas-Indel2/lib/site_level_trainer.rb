# #!/data/pmrs/jin/bin/ruby -I /data/pmrs/challis/lib/Atlas-Indel2/trunk/
$:.unshift File.join(File.dirname(__FILE__))
require 'read.rb'
require 'read_indel.rb'
require 'read_insertion.rb'
require 'read_deletion.rb'
require 'assertion_failure.rb'

# #require 'getoptlong'
require 'getoptlong.rb'
require 'ref_seq.rb'
require 'indel.rb'

$DEBUG=false
@@e=2.71828182845904523536
@@buffer_size = 300000
@@illum_read_level_z_cutoff = -4.7 # illumina read cutoff
# put read level model coefficients here:
@@illum_read_intercept = -8.17632
#@@illum_read_average_nbq = 0.08582
@@illum_read_var_rate = -74.97455
@@illum_read_mean_base_quality = 0.21481
@@illum_read_is_clipped = 0.77915

@@solid_read_level_z_cutoff = -2.5 # SOLiD read cutoff
@@solid_read_intercept = -6.30828
@@solid_read_average_nbq = 0.13705
@@solid_read_var_rate = -1.90628
#@@solid_read_near_end = -0.929674

class Site_level_trainer
	def initialize(infile, ref_filename, insertion_outfile, deletion_outfile, illumina, heuristics=false, true_positives=true)
		if(true_positives)
			@true_indel = 1
		else
			@true_indel = 0
		end
		@heur=heuristics
		if(illumina)
			@platform = "illumina"
			@read_level_z_cutoff = @@illum_read_level_z_cutoff
		else
			@platform = "solid"
			@read_level_z_cutoff = @@solid_read_level_z_cutoff
		end
		begin
			@reference = Ref_seq.new(ref_filename)
		rescue
			raise "Error loading the reference sequence from file #{ref_filename}:\n#{$!}"
		end
		begin
			@num_lines = `wc -l #{infile}`.to_i
			@infile = File.open(infile, 'r')
			@insertion_outfile = File.open(insertion_outfile, 'w')
			@deletion_outfile = File.open(deletion_outfile, 'w')
		rescue
			raise "There was an error loading the input or output file:\n#{$!}"
		end
	end

	def self.process_arguments()
		opts = GetoptLong.new(
		["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
		["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT ],
		["--insertion-outfile","-I", GetoptLong::REQUIRED_ARGUMENT],
		["--deletion-outfile","-D", GetoptLong::REQUIRED_ARGUMENT],
		["--true-negatives","-n", GetoptLong::NO_ARGUMENT],
		["--heuristics", "-h", GetoptLong::NO_ARGUMENT],
		["-S", GetoptLong::NO_ARGUMENT]

		)

		opt_hash = {}
		opts.each do |opt, arg|
			opt_hash[opt] = arg
		end

		if( !(opt_hash.key?("--input")) || !(opt_hash.key?("--deletion-outfile")) || !(opt_hash.key?("--reference")) || !(opt_hash.key?("--insertion-outfile")) )
			puts "USAGE:\n	site_level_trainer.rb -i infile -r reference_sequence -I insertion_outfile -D deletion_outfile --true-negatives -h (heuristics) -S (SOLiD)\n"
			exit 0
		end

		return opt_hash
	end

	def process_data
		skipped = 0
		skipped_flagged = 0
		processed = 0
		progress = 0
		buffered_coor = 0 #how far has the current chromosome been buffered?
		@insertion_outfile.puts("true_indel\tgeno\tchr\tcoor\tref\talt\talt_depth\ttotal_depth\tpassed_reads\tvar_read_ratio\tpassed_read_ratio\tmean_z\tstrand_dir\tframeshift\tlength\tmean_map_qual\tmean_ave_nqs\tmean_var_rate\tlocal_entropy\tsimple_local_entropy\tnorm_var_square\tsample\tstrand_score\tnew_strand_filt\tnear_read_end\tread_end_ratio\tread_end_plus_strand\tmean_color_corrections\tmax_map_qual\thigh_var_ratio\tmean_var_rate2\tartifact_variant_ratio\tartifact_variants")
		@deletion_outfile.puts("true_indel\tgeno\tchr\tcoor\tref\talt\talt_depth\ttotal_depth\tpassed_reads\tvar_read_ratio\tpassed_read_ratio\tmean_z\tstrand_dir\tframeshift\tlength\tmean_map_qual\tmean_ave_nqs\tmean_var_rate\tlocal_entropy\tsimple_local_entropy\tnorm_var_square\tsample\tstrand_score\tnew_strand_filt\tnear_read_end\tread_end_ratio\tread_end_plus_strand\tmean_color_corrections\tmax_map_qual\thigh_var_ratio\tmean_var_rate2\tartifact_variant_ratio\tartifact_variants")
		chr = nil
		indel_buffer = Hash.new
		depth = Hash.new
		line = get_line()
		geno = './.'
		while(line != :eof)
			begin
				if(line =~ /^###/)
					# end of this indel, print out and reset and move on
					print_buffer(indel_buffer, 1.0/0.0, chr, depth) #pass infinity as the buffered_coor to empty the buffer
					indel_buffer = Hash.new
					depth = Hash.new
					next
				end
				if(line =~ /SAMP:Z:([^	]+)/)
					sample = $1
				end
				read = Read.new(line, false)
				chr = read.ref if chr == nil
				buffered_coor = read.pos
				if(indel_buffer.length > @@buffer_size && chr == read.ref) #buffer is full, print buffer
					print_buffer(indel_buffer, buffered_coor, chr, depth)
				elsif(chr != read.ref) #reached end of chromosome, print buffer
					print_buffer(indel_buffer, 1.0/0.0, chr, depth) #pass infinity as the buffered_coor to empty the buffer
					assert(indel_buffer.length == 0, "ERROR: Indel buffer is not empty after printing end of chromosome #{chr}. Buffer length: #{indel_buffer.length}\n#{indel_buffer.to_s}")
					indel_buffer = Hash.new
					chr = read.ref
					depth = Hash.new
				end
				# update depth info
				indels = read.get_indels # this must be run before read.end_pos can be called
				(read.pos..read.end_pos).each do |coor|
					if(depth[coor] == nil)
						depth[coor] = 1
					else
						depth[coor]+= 1
					end
				end
				if(line =~ /COOR:i:(\d+)/)
					tp_coor = $1.to_i
				else
					tp_coor = 0
				end
				if(line =~ /GENO:Z:(...)/)
					geno = $1
				end
				indels.each do |read_indel|
					c = tp_coor
					if(c == read_indel.start_pos) # make sure it is the real TP indel
						process_indel(read_indel, indel_buffer, read, sample, geno)
					end
				end
				processed += 1
			rescue FormatError
				puts $!.to_s + "\n\n" + line + "\n\n"
				skipped += 1
			rescue FlaggedError
				skipped_flagged += 1
				skipped += 1
				# #rescue # puts "There was an unkown error with a line, it will be skipped
				# \n\n#{line}\n\n#{$!}" # skipped += 1
			ensure
				if( progress != ((processed+skipped).to_f/@num_lines.to_f * 100).to_i )
					progress = indicate_progress(processed, skipped)
				end
				line = get_line()
			end
		end
		@infile.close
		print "writing output from buffer..."
		STDOUT.flush
		print_buffer(indel_buffer, 1.0/0.0, chr, depth)  #print whatever is left in the buffer
		assert(indel_buffer.length == 0, "ERROR: Indel buffer is not empty at the end of the last chromosome. Buffer length: #{indel_buffer.length}\n#{indel_buffer.to_s}")
		STDOUT.flush
		@insertion_outfile.close
		@deletion_outfile.close
		puts "#{skipped} lines skipped"
		puts "#{skipped_flagged} flagged lines skipped"
		puts "#{processed} lines processed"
	end

	private

	def indicate_progress(processed, skipped)
		puts "#{((processed+skipped).to_f/@num_lines.to_f * 100).to_i}%"
		STDOUT.flush
		STDERR.flush
		@deletion_outfile.flush
		@insertion_outfile.flush
		return ((processed+skipped).to_f/@num_lines.to_f * 100).to_i
	end

	def is_minor_allele(buffer, indel_key)
		indel = buffer[indel_key][0]
		buffer.each do |key, value|
			next if key == indel_key
			if(indel.read_indel.start_pos == value[0].read_indel.start_pos && indel.read_count <= value[0].read_count)
				return true
			end
		end
		return false
	end

	def process_indel(read_indel, indel_buffer, read, sample, geno)
		z=nil
		z = read_level_model(read_indel.avnbq, read_indel.near_read_end, read_indel.read.mean_base_qual, read_indel.read.is_clipped)
		indel_key = read_indel.to_s
		if( indel_buffer[indel_key] == nil)
			indel_buffer[indel_key] = [Indel.new(read_indel, z, @read_level_z_cutoff, sample), geno]
		else
			indel_buffer[indel_key][0].add_read_indel(z, read_indel)
		end
	end

	def print_buffer(buffer, buffer_coor, chr, depth)
		print "writing output from buffer..."
		STDOUT.flush
		i=buffer.length
		buffer.delete_if {|key,value| is_minor_allele(buffer, key)} # remove minor alleles
		buffer.each do |key, value|
			indel = value[0]
			geno = value[1]
			if(indel.read_indel.start_pos < buffer_coor)

				raise "ERROR: missing total depth information for indel #{indel.read_indel.to_s} -- #{indel.read_indel.start_pos()}" if depth[indel.read_indel.start_pos()] == nil
				assert( chr == key.split(":")[0].to_sym, "ERROR: tried to print #{indel} out on chromosome #{chr}")
				assert( key.split(/[:DI]/)[1].to_i == indel.read_indel.start_pos, "ERROR: tried to print #{indel} out at wrong coordinate")
				buffer.delete(key)
				var_reads_ratio = ((indel.read_count.to_f / depth[indel.read_indel.start_pos()].to_f)*1000).round/1000.0
				passed_reads_ratio = ((indel.num_passed_reads.to_f / depth[indel.read_indel.start_pos()].to_f)*1000).round/1000.0
				if(indel.fails_strand_dir_filt)
					strand_dir = 0
				else
					strand_dir= 1
				end
				if(indel.read_indel.is_deletion)
					# #@insertion_outfile.puts("chr\tcoor\tref\talt\talt_depth\ttotal_depth\tpassed_reads\tvar_read_ratio\tpassed_read_ratio\tmean_z\tsum_z\tframeshift\tlength\tmean_map_qual\tmean_ave_nqs\tmean_var_rate\tlocal_entropy\tmean_simple_local_entropy")
					@deletion_outfile.puts "#{@true_indel}\t#{geno}\t#{chr}\t#{indel.read_indel.start_pos}\t#{indel.read_indel.ref_seq}\t#{indel.read_indel.alt_seq}\t#{indel.read_count}\t#{depth[indel.read_indel.start_pos()]}\t#{indel.num_passed_reads}\t#{var_reads_ratio}\t#{passed_reads_ratio}\t#{indel.z_score}\t#{strand_dir}\t#{indel.frameshift}\t#{indel.length}\t#{indel.mean_map_qual}\t#{indel.mean_avnqs}\t#{indel.mean_var_rate}\t#{indel.local_entropy}\t#{indel.simple_local_entropy}\t#{(indel.read_count**2).to_f/depth[indel.read_indel.start_pos()].to_f}\t#{indel.sample}\t#{indel.strand_score}\t#{indel.new_strand_dir}\t#{indel.all_near_read_end}\t#{indel.near_read_end_ratio}\t#{indel.near_read_end_and_strand_bais}\t#{indel.mean_color_corrections}\t#{indel.max_map_qual}\t#{indel.high_var_ratio}\t#{indel.mean_var_rate2}\t#{indel.artifact_variant_ratio}\t#{indel.artifact_variants}"
				else
					@insertion_outfile.puts "#{@true_indel}\t#{geno}\t#{chr}\t#{indel.read_indel.start_pos}\t#{indel.read_indel.ref_seq}\t#{indel.read_indel.alt_seq}\t#{indel.read_count}\t#{depth[indel.read_indel.start_pos()]}\t#{indel.num_passed_reads}\t#{var_reads_ratio}\t#{passed_reads_ratio}\t#{indel.z_score}\t#{strand_dir}\t#{indel.frameshift}\t#{indel.length}\t#{indel.mean_map_qual}\t#{indel.mean_avnqs}\t#{indel.mean_var_rate}\t#{indel.local_entropy}\t#{indel.simple_local_entropy}\t#{(indel.read_count**2).to_f/depth[indel.read_indel.start_pos()].to_f}\t#{indel.sample}\t#{indel.strand_score}\t#{indel.new_strand_dir}\t#{indel.all_near_read_end}\t#{indel.near_read_end_ratio}\t#{indel.near_read_end_and_strand_bais}\t#{indel.mean_color_corrections}\t#{indel.max_map_qual}\t#{indel.high_var_ratio}\t#{indel.mean_var_rate2}\t#{indel.artifact_variant_ratio}\t#{indel.artifact_variants}"
				end
				@deletion_outfile.flush
				@insertion_outfile.flush
			end
		end
		depth.delete_if { |coor, depth| coor < buffer_coor } #clean out the depth buffer
		print "#{i-buffer.length} lines printed\n"
		STDOUT.flush
	end

	def read_level_model(average_nbq, near_read_end, mean_base_quality, is_clipped)
		return 0
		#		if(is_clipped)
		#			is_clipped=1.0
		#		else
		#			is_clipped=0.0
		#		end
		#		case @platform
		#			when "illumina"
		#				return @@illum_read_intercept + @@illum_read_mean_base_quality*mean_base_quality.to_f + @@illum_read_var_rate*var_rate + @@illum_read_is_clipped*is_clipped
		#			when "solid"
		#				return @@solid_read_intercept + @@solid_read_average_nbq*average_nbq.to_f + @@solid_read_var_rate*var_rate
		#			else
		#				raise "missing or unknown sequencing platform"
		#		end

	end

	def get_line
		begin
			line = @infile.readline
			line.chomp!
		rescue EOFError => err
			return :eof
		end
		return line
	end

end

# main
if __FILE__ == $0
	opt_hash = Site_level_trainer.process_arguments()
	atlas = Site_level_trainer.new(opt_hash["--input"], opt_hash["--reference"], opt_hash["--insertion-outfile"],  opt_hash["--deletion-outfile"], opt_hash["-S"]==nil, opt_hash["--heuristics"]!=nil,  opt_hash["--true-negatives"]==nil)
	atlas.process_data

	print "DONE\n\nAtlas-INDEL2 site-level trainer Finished\n"
end
