# #this class is responsible for pulling true positives, as defined in a VCF
# input file, (one sample-site per line) #out of a list of bam files to create a
# file of true-positive indel SAM lines #to each sam line is appened the fields
# CHR:Z:[chromosome] and COOR:i:[coordinate] to specify #the exact location of
# the indel in the read #BCM-HGSC #Author: Danny Challis, Fuli Yu #30 Aug 2010
$:.unshift File.join(File.dirname(__FILE__))
class True_positive_site_generator
	require 'data_file'
	require 'read'

	def initialize( tp_file, bam_file_array, allow_multi_sample, true_positives)
		@allow_multi_sample = allow_multi_sample
		@tp_file = Data_file.new(tp_file) #true positives file
		@bam_files = bam_file_array #source bam files
		@true_positives = true_positives # is this true positive data or false positive data?
		# #@bam_file = bam_file
		STDERR.puts "Allowing multiple samples per loci" if @allow_multi_sample
	end


	def run(writer)
		@tp_file.process do |tp_cols|  #reads through file, splits into array of data columns (true-positive_columns)
			if(tp_cols[0][0,6]=='#CHROM')#does NOT check for site filters
				@header=tp_cols
				next
			elsif tp_cols[0][0,1]=='#'
				next
			end
			chr = tp_cols[0]
			coor = tp_cols[1].to_i
			ref = tp_cols[3]
			alts = tp_cols[4].split(",")
			#length = tp_cols[4].length-1 
			if(@allow_multi_sample)
				samples,genotypes = get_samples(tp_cols)
			else
				sample, geno = get_sample(tp_cols)    # these lines
				samples = [sample]
				genotypes = [geno]
			end
			if(samples.length == 0)
				STDERR.puts "No samples found for #{chr}:#{coor}"
			else 
				STDERR.puts "Found #{samples.length} samples for #{chr}:#{coor}" if $DEBUG
			end
			coor = coor+1 if tp_cols[3].length > 1
			samples = samples.sort_by {rand} #shuffle the order of the samples list
			samples.each_with_index do |sample, i|
				geno = genotypes[i]
				if(geno =~ /[1-9]/)
					allele_i = $&.to_i - 1
				end
				alt = alts[allele_i]
				selected_lines = Hash.new
				output = Array.new
				ref_lines=Array.new
				bam_file = get_bam_file(sample, chr)
				if bam_file == nil # file is missing, skip it
					STDERR.puts "could not find BAM file for #{sample}, skipping..."
					next
				end
				sam_lines_str = `samtools view -F 1796 #{bam_file} #{chr}:#{coor-50}-#{coor+50}`
				sam_lines = sam_lines_str.split("\n")
				# sam_lines = sam_lines.sort_by {rand} #shuffle the order of the sam lines
				pileup = Hash.new
				sam_lines.each do |read|
					begin
						nearby_coor, a_length = read_has_indel(read, tp_cols)
						if(nearby_coor != nil)
							if(pileup["#{nearby_coor}_#{a_length}"] == nil)
								pileup["#{nearby_coor}_#{a_length}"] = Array.new
							end
							pileup["#{nearby_coor}_#{a_length}"].push( "#{read}	SAMP:Z:#{sample}	ALLELE:Z:#{indel_key}	GENO:Z:#{geno}" )
							nearby_coor = nil
						else
							ref_lines.push("#{read}	SAMP:Z:#{sample}")
						end
					rescue FlaggedError
						# skip this line
					end
				end
				max = 0
				consensus = nil
				consensus_key = nil
				pileup.each do |key, sam_arr|
					if sam_arr.length > max
						max = sam_arr.length
						consensus = sam_arr
						consensus_key = key
					end
				end
				if @true_positives && max < 2 # require at least 2 variant reads for TP data
					selected_lines = Hash.new
					output = Array.new
					next
				end
				if(consensus != nil ) 
					consensus.each do |line|
						output.push(line)
						selected_lines[consensus_key] = 0 if selected_lines[consensus_key].nil?
						selected_lines[consensus_key] += 1
					end
					ref_lines.each do |line| # include reference reads
						output.push(line)
					end
					pileup.delete(consensus_key) # make sure to include any nearby gapped reads that didn't make the cut
					pileup.each do |key, sam_arr|
						sam_arr.each do |line|
							if(line =~ /(.*)	SAMP:Z:/)
								output.push($1)
							end
						end
					end
					if(selected_lines == 0)
						STDERR.puts "Indel at #{chr}:#{coor} not found in any #{sample}"
					else
						output.sort! {|a,b| a.split("\t")[3].to_i <=> b.split("\t")[3].to_i }
						output.each {|outline| writer.puts outline}
						writer.puts '###########################################' unless output.length == 0
					end
					break if( !@allow_multi_sample && selected_lines[consensus_key] >= 1 ) # stop after you find at least on line for a sample (we don't want multiple samples per site)
				end	# if(consensus
			end # samples.each
		end # @tp_file.process do
	end # def


	private #--------------------------------------------------------------

	# get the list of samples from the realignment interval files -- this is
	# specific to pilot 3 training
#	def get_samples(chr, coor)
#		samples = Array.new
#		grep_result = `grep "#{chr}:#{coor-100}" intervals/*.output.intervals`
#		grep_result.split("\n").each do |interval_file|
#			samples.push(interval_file.slice(/NA\d+/))
#		end
#		return samples
#	end

	# verify that the read has the validated indel, allows for a 5bp discrepency
	def read_has_indel(read_str, cols)
		coor = cols[1].to_i
		if(cols[3].length == 1)
			event = :"I"
			 length = cols[4].length - 1
		else
			event = :"D"
			 length = cols[3].length - 1
		end
		read = Read.new(read_str, false)
		index = read.pos - 1
		read.cigar.each do |entry|
			if( (index - coor).abs < 10 )
				if(entry[1] == event   && (entry[0]-length).abs < 1)
					return index, entry[0]
				end
			end
			index += entry[0] unless entry[1] == :"I"
		end
		return nil
	end

	# #get the validated sample names -- LIMITATION: does not check for individual
	# sample filters
	def get_sample(cols)
		if(cols[7] =~ /sample=([^;\t]+)/)#search info column for validated sample name
			sample = $1
			if(cols[9] =~ /([1234]\/.)|(.\/[1234])/)
				return sample, $&
			else
				return sample, "./."
			end
		end
		if(@header == nil)
			raise "Cannot determine sample name, no \"sample\" entry in info column and no VCF header provided."
		end
		index_array = (9...(cols.length))
		index_array = index_array.sort_by {rand} #shuffle the samples to be more random
		index_array.each do |i| #search sample columns for validated sample
			if(cols[i] =~ /([1234]\/.)|(.\/[1234])/)
				return @header[i], $&
			end
		end
	end

	# #get the variant sample names -- LIMITATION: does not check for individual
	# sample filters
	def get_samples(cols)
		samples = Array.new
		genotypes = Array.new
		if(@header == nil)
			raise "Cannot determine sample name, no VCF header provided."
		end
		index_array = (9...(cols.length))
		index_array = index_array.sort_by {rand} #shuffle the samples to be more random
		index_array.each do |i| #search sample columns for validated sample
			if(cols[i] =~ /([1234]\/.)|(.\/[1234])/)
				samples.push @header[i]
				genotypes.push($&)
			end
		end
		return samples,genotypes
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


end

if(ARGV.length < 2)
	puts "Usage: true_positive_generator.rb [tp_vcf_file] \"path/to/bam/files/*.bam\" T/F (true positive or false positive data?) -m (allow multiple samples per loci)"
	exit 0
end

tp_file = ARGV[0]
bam_files = Dir.glob(ARGV[1])
# bam_files = ARGV[1]
STDERR.puts "Starting up..."
tp_gen = True_positive_site_generator.new(tp_file, bam_files, ARGV[3]=='-m', ARGV[2]=='T')
STDERR.puts "running..."
tp_gen.run(STDOUT)
STDERR.puts "DONE"
