# #!/data/pmrs/jin/bin/ruby -I /data/pmrs/challis/lib/Atlas-Indel2/trunk/
$:.unshift File.join(File.dirname(__FILE__))
require 'read.rb'
require 'vcf_line.rb'
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

class Atlas_indel_trainer
	def initialize(bam_files, vcf_file, ref_filename, insertion_outfile, deletion_outfile, true_positives_only, read_level)
		@true_only = true_positives_only
		@read_level = read_level
		begin
			@reference = Ref_seq.new(ref_filename)
		rescue
			raise "Error loading the reference sequence from file #{ref_filename}:\n#{$!}"
		end
		begin
			@bam_files = bam_files
			@vcf_file = vcf_file
			@insertion_outfile = File.open(insertion_outfile, 'w')
			@deletion_outfile = File.open(deletion_outfile, 'w')
		rescue
			raise "There was an error loading the input or output file:\n#{$!}"
		end
	end

	def self.process_arguments()
		opts = GetoptLong.new(
		["--bams", "-b", GetoptLong::REQUIRED_ARGUMENT],
		["--vcf", "-v", GetoptLong::REQUIRED_ARGUMENT],
		["--reference", "-r", GetoptLong::REQUIRED_ARGUMENT ],
		["--insertion-outfile","-I", GetoptLong::REQUIRED_ARGUMENT],
		["--deletion-outfile","-D", GetoptLong::REQUIRED_ARGUMENT],
		["--true-negatives","-n", GetoptLong::NO_ARGUMENT],
		["--read-level","-R", GetoptLong::NO_ARGUMENT]

		)

		opt_hash = {}
		opts.each do |opt, arg|
			opt_hash[opt] = arg
		end

		if( !(opt_hash.key?("--bams")) || !(opt_hash.key?("--vcf")) || !(opt_hash.key?("--deletion-outfile")) || !(opt_hash.key?("--reference")) || !(opt_hash.key?("--insertion-outfile")) )
			puts "USAGE:\n	atlas_indel_trainer.rb -b \"*.bam\" -v file.vcf -r reference_sequence -I insertion_outfile -D deletion_outfile --true-negatives --read-level\n"
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
		if(@read_level)
			@insertion_outfile.puts("true_indel\tgeno\tchr\tcoor\tref\talt\tframeshift\tlength\tmap_qual\tnqs\tvar_rate\tlocal_entropy\tsimple_local_entropy\tsample\tnear_read_end\tcolor_corrections\tvar_count\tmean_base_qual\tdist_3")
			@deletion_outfile.puts("true_indel\tgeno\tchr\tcoor\tref\talt\tframeshift\tlength\tmap_qual\tnqs\tvar_rate\tlocal_entropy\tsimple_local_entropy\tsample\tnear_read_end\tcolor_corrections\tvar_count\tmean_base_qual\tdist_3")
		else
			@insertion_outfile.puts("true_indel\tgeno\tchr\tcoor\tref\talt\talt_depth\ttotal_depth\tref_depth\tvar_read_ratio\tstrand_dir\tframeshift\tlength\tmean_map_qual\tmean_ave_nqs\tmean_var_rate\tlocal_entropy\tsimple_local_entropy\tnorm_var_square\tsample\tstrand_score\tnew_strand_filt\tnear_read_end\tread_end_ratio\tread_end_plus_strand\tmean_color_corrections\tmax_map_qual\thigh_var_ratio\tmean_var_rate2\tartifact_variant_ratio\tartifact_variants")
			@deletion_outfile.puts("true_indel\tgeno\tchr\tcoor\tref\talt\talt_depth\ttotal_depth\tref_depth\tvar_read_ratio\tstrand_dir\tframeshift\tlength\tmean_map_qual\tmean_ave_nqs\tmean_var_rate\tlocal_entropy\tsimple_local_entropy\tnorm_var_square\tsample\tstrand_score\tnew_strand_filt\tnear_read_end\tread_end_ratio\tread_end_plus_strand\tmean_color_corrections\tmax_map_qual\thigh_var_ratio\tmean_var_rate2\tartifact_variant_ratio\tartifact_variants")
		end
		chr = nil

		labels=nil
		File.open(@vcf_file, 'r').each_line do |vcf_str|
			vcf_str.chomp!
			if(vcf_str[0,1] == '#')
				if( vcf_str[0,4] == '#CHR' )
					labels = vcf_str.split("\t")
				end
				next
			end
			vcf_line = Vcf_line.read_line(vcf_str, false, labels)
			ref = vcf_line.ref().to_s
			vcf_line.sample_names.each do |sample_name|
				geno = vcf_line.samples[sample_name][:"GT"].to_s
				if(geno =~ /[1-9]/)
					allele_i = $&.to_i - 1
				elsif( geno =~ /0/ && !@true_only)
					allele_i = 0 # limitation, only checks for the first allele true negative
				else
					next
				end
				alt = vcf_line.alt()[allele_i].to_s
				if( alt.length > ref.length )
					ins_bases = alt[1..(alt.length-ref.length)]
					indel_key = "#{vcf_line.chr}:#{vcf_line.pos}I#{ins_bases}"
				else
					indel_key = "#{vcf_line.chr}:#{vcf_line.pos}D#{ref.length - alt.length}"
				end

				bam_file =get_bam_file(sample_name, vcf_line.chr())
				if(bam_file.nil?)
					STDERR.puts "skipping to next sample..."
					next
				end
				@infile = IO.popen("samtools view -F 1796 #{bam_file} #{vcf_line.chr}:#{vcf_line.pos-1}-#{vcf_line.pos+1}")
				line = get_line()
				indel_buffer = Hash.new
				depth = Hash.new
				ref_depth = Hash.new
				the_indel = nil # this will hold the indel object we are searching for
				while(line != :eof)
					begin
						read = Read.new(line, false)
						chr = read.ref

						# update depth info
						indels = read.get_indels_with_depth(depth, ref_depth)

						# don't be picky about the allele for true negatives
						if(!@true_only && geno.to_s == "0/0")
							indels.each do |read_indel|
								if(read_indel.start_pos == vcf_line.pos)
									indel_key = read_indel.to_s
									break
								end
							end
						end

						# gather indel reads
						indels.each do |read_indel|
							if(read_indel.to_s == indel_key) # make sure it's the right indel
								if(@read_level)
									print_read(read_indel, sample_name, geno)
									the_indel = :found_it
								else
									if(the_indel.nil?)
										the_indel = Indel.new(read_indel, 0, 0, sample_name)
									else
										the_indel.add_read_indel(0,read_indel)
									end
								end
							end
						end
					rescue FormatError
						puts $!.to_s + "\n\n" + line + "\n\n"
						skipped += 1
					rescue FlaggedError
						skipped_flagged += 1
						skipped += 1
						# #rescue # puts "There was an unkown error with a line, it will be skipped
						# \n\n#{line}\n\n#{$!}" # skipped += 1
					ensure
						line = get_line()
					end
				end
				@infile.close
				if the_indel.nil?
					STDERR.puts "Could not find #{indel_key} in sample #{sample_name}"
					STDERR.flush
					next
				end
				print_indel(the_indel, vcf_line.chr(), geno, depth, ref_depth) unless @read_level
			end
		end
		@insertion_outfile.close
		@deletion_outfile.close
	end

	private

	def print_read(read_indel, sample, geno)
		if(geno.to_s == "0/0")
			is_true = 0
		else
			is_true = 1
		end
		if(read_indel.is_deletion)
			@deletion_outfile.puts("#{is_true}\t#{geno}\t#{read_indel.chr}\t#{read_indel.start_pos}\t#{read_indel.ref_seq}\t#{read_indel.alt_seq}\t#{read_indel.frame_shift}\t#{read_indel.length}\t#{read_indel.map_qual}\t#{read_indel.avnbq}\t#{read_indel.var_rate}\t#{read_indel.local_entropy}\t#{read_indel.simple_local_entropy}\t#{sample}\t#{read_indel.near_read_end}\t#{read_indel.read.color_correction}\t#{read_indel.read.var_count}\t#{read_indel.read.mean_base_qual}\t#{read_indel.dist_3}")
		else
			@insertion_outfile.puts("#{is_true}\t#{geno}\t#{read_indel.chr}\t#{read_indel.start_pos}\t#{read_indel.ref_seq}\t#{read_indel.alt_seq}\t#{read_indel.frame_shift}\t#{read_indel.length}\t#{read_indel.map_qual}\t#{read_indel.avnbq}\t#{read_indel.var_rate}\t#{read_indel.local_entropy}\t#{read_indel.simple_local_entropy}\t#{sample}\t#{read_indel.near_read_end}\t#{read_indel.read.color_correction}\t#{read_indel.read.var_count}\t#{read_indel.read.mean_base_qual}\t#{read_indel.dist_3}")
		end
	end

	# #searches the @bam_files array for a bam file with the needed sample name
	def get_bam_file(sample, chr)
		@bam_files.each do |filename|
			if(filename.include?(sample))
				if(filename.include?('chrom'))
					if(filename.include?("chrom#{chr}."))
						return filename
					end
				else
					return filename
				end
			end
		end
		STDERR.puts "Could not find BAM file for #{sample} at chromosome #{chr}"
		return nil
	end

	def indicate_progress(processed, skipped)
		puts "#{((processed+skipped).to_f/@num_lines.to_f * 100).to_i}%"
		STDOUT.flush
		STDERR.flush
		@deletion_outfile.flush
		@insertion_outfile.flush
		return ((processed+skipped).to_f/@num_lines.to_f * 100).to_i
	end

	def print_indel(indel, chr, geno, depth, ref_depth)

		raise "ERROR: missing total depth information for indel #{indel.read_indel.to_s} -- #{indel.read_indel.start_pos()}" if depth[indel.read_indel.start_pos()] == nil
		raise "ERROR: missing ref depth information for indel #{indel.read_indel.to_s} -- #{indel.read_indel.start_pos()}" if ref_depth[indel.read_indel.start_pos()] == nil
		var_reads_ratio = ((indel.read_count.to_f / depth[indel.read_indel.start_pos()].to_f)*1000).round/1000.0
		passed_reads_ratio = ((indel.num_passed_reads.to_f / depth[indel.read_indel.start_pos()].to_f)*1000).round/1000.0
		if(indel.fails_strand_dir_filt)
			strand_dir = 0
		else
			strand_dir= 1
		end
		if(geno.to_s == "0/0")
			is_true = 0
		else
			is_true = 1
		end
		raise "missing genotype" if geno.to_s == "" || geno.to_s == "./."
		if(indel.read_indel.is_deletion)
			# #@insertion_outfile.puts("chr\tcoor\tref\talt\talt_depth\ttotal_depth\tpassed_reads\tvar_read_ratio\tpassed_read_ratio\tmean_z\tsum_z\tframeshift\tlength\tmean_map_qual\tmean_ave_nqs\tmean_var_rate\tlocal_entropy\tmean_simple_local_entropy")
			@deletion_outfile.puts "#{is_true}\t#{geno}\t#{chr}\t#{indel.read_indel.start_pos}\t#{indel.read_indel.ref_seq}\t#{indel.read_indel.alt_seq}\t#{indel.read_count}\t#{depth[indel.read_indel.start_pos()]}\t#{ref_depth[indel.read_indel.start_pos()]}\t#{var_reads_ratio}\t#{strand_dir}\t#{indel.frameshift}\t#{indel.length}\t#{indel.mean_map_qual}\t#{indel.mean_avnqs}\t#{indel.mean_var_rate}\t#{indel.local_entropy}\t#{indel.simple_local_entropy}\t#{(indel.read_count**2).to_f/depth[indel.read_indel.start_pos()].to_f}\t#{indel.sample}\t#{indel.strand_score}\t#{indel.new_strand_dir}\t#{indel.all_near_read_end}\t#{indel.near_read_end_ratio}\t#{indel.near_read_end_and_strand_bais}\t#{indel.mean_color_corrections}\t#{indel.max_map_qual}\t#{indel.high_var_ratio}\t#{indel.mean_var_rate2}\t#{indel.artifact_variant_ratio}\t#{indel.artifact_variants}"
		else
			@insertion_outfile.puts "#{is_true}\t#{geno}\t#{chr}\t#{indel.read_indel.start_pos}\t#{indel.read_indel.ref_seq}\t#{indel.read_indel.alt_seq}\t#{indel.read_count}\t#{depth[indel.read_indel.start_pos()]}\t#{ref_depth[indel.read_indel.start_pos()]}\t#{var_reads_ratio}\t#{strand_dir}\t#{indel.frameshift}\t#{indel.length}\t#{indel.mean_map_qual}\t#{indel.mean_avnqs}\t#{indel.mean_var_rate}\t#{indel.local_entropy}\t#{indel.simple_local_entropy}\t#{(indel.read_count**2).to_f/depth[indel.read_indel.start_pos()].to_f}\t#{indel.sample}\t#{indel.strand_score}\t#{indel.new_strand_dir}\t#{indel.all_near_read_end}\t#{indel.near_read_end_ratio}\t#{indel.near_read_end_and_strand_bais}\t#{indel.mean_color_corrections}\t#{indel.max_map_qual}\t#{indel.high_var_ratio}\t#{indel.mean_var_rate2}\t#{indel.artifact_variant_ratio}\t#{indel.artifact_variants}"
		end
		@deletion_outfile.flush
		@insertion_outfile.flush
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
	opt_hash = Atlas_indel_trainer.process_arguments()
	atlas = Atlas_indel_trainer.new(Dir.glob(opt_hash["--bams"]), opt_hash["--vcf"], opt_hash["--reference"], opt_hash["--insertion-outfile"],  opt_hash["--deletion-outfile"], opt_hash["--true-negatives"]==nil, opt_hash["--read-level"]!=nil)
	atlas.process_data

	print "DONE\n\natlas_indel_trainer Finished\n"
end
