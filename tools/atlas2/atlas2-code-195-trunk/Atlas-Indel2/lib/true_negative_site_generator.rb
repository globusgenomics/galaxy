$:.unshift File.join(File.dirname(__FILE__))
# CHR:Z:[chromosome] and COOR:i:[coordinate] to specify #the exact location of
# the indel in the read #BCM-HGSC #Author: Danny Challis, Fuli Yu #30 Aug 2010
class True_negative_site_generator
	require 'data_file'
	require 'read'

	def initialize( pileup_file_array, bam_file_array, call_set, solid_platform)
		@pileup_files = pileup_file_array
		@bam_files = bam_file_array #source bam files
		@call_set = call_set
		@solid = solid_platform
	end


	def run()
		@included_sites = Hash.new
		@pileup_files.each do |pileup_file|
			if(pileup_file =~ /\/([^\.\/]+)[^\/]*$/)
				sample = $1
			else
				raise "no sample name in pileup file name!  \n #{pileup_file}"
			end
			found_lines = 0
			File.open(pileup_file,'r').each_line do |line|
				break if found_lines > 50
				next if rand > 0.001 # distribute across the genome
				@included_reads = Hash.new
				next if( rand(3) != 1) # take a random selection, not all the data as there's too much
				line.chomp!
				next if line == "" || line == nil
				if(line =~ /^([\dXY]+)	(\d+)/)
					chr = $1
					coor = $2.to_i
				else
					raise "malformed pileup line: \n #{line}"
				end
				next if !@included_sites["#{chr}:#{coor}"].nil? # this site was already included
				@included_sites["#{chr}:#{coor}"] = sample # this site was already included
				if(chr =~ /\d+/)
					chr_code = (chr.to_i-1).to_s(27).tr("0-9a-q", "A-Z") # make alphabetical chr code
				else
					chr_code = chr
				end
				writer = File.open("tn_sites.#{chr_code}", 'a')
				# writer = File.open("tn_sites.sam", 'a') # don't need separate chromosome files this time
				if( `awk '(($2-#{coor})^2)^0.5 <= 5' #{@call_set}  | wc -l`.to_i > 0 ) # skip lines in the call set
					STDERR.puts "Indel at #{chr}:#{coor} #{sample} is in the call set, skipping..."
					next
				end
				bam_file = get_bam_file(sample, chr)
				next if bam_file == nil # file is missing, skip it
				sam_lines_str = `samtools view #{bam_file} #{chr}:#{coor-50}-#{coor+50}`
				sam_lines = sam_lines_str.split("\n")
				selected_lines = 0
				output = Array.new
				sam_lines.each do |read|
					begin
						nearby_coor, a_length = read_has_indel(read, coor, sample)
						if(nearby_coor != nil) # read has the indel
							output << "#{read}	SAMP:Z:#{sample}	CHR:Z:#{chr}	COOR:i:#{nearby_coor}\n"
							selected_lines += 1
							nearby_coor = nil
						else #read does not have the indel
							output << "#{read}	SAMP:Z:#{sample}\n"
						end
					rescue FlaggedError
						# skip line
					end
				end
				if(selected_lines == 0)
					STDERR.puts "Indel at #{chr}:#{coor} not found in #{sample}"
				else
					output.each {|outline| writer.puts outline}
					writer.puts '###########################################'
					found_lines += 1
				end
				writer.close
			end
		end
	end


	private #--------------------------------------------------------------


	# #verify that the read has the validated indel, allows for a 5bp discrepency
	def read_has_indel(read_str, coor, sample)
		read = Read.new(read_str, @solid)
		raise FlaggedError if !@included_reads["#{read.name}_#{read.flag}"].nil? # this read was already included
		@included_reads["#{read.name}_#{read.flag}"]=sample
		index = read.pos - 1
		read.cigar.each do |entry|
			if( index == coor && (entry[1] == :"I" || entry[1] == :"D") )
				return index, entry[0]
			end
			index += entry[0] unless entry[1] == :"I"
		end
		return nil
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

if(ARGV.length < 3)
	puts "Usage: true_negative_generator.rb \"path/to/pileup/files/*.pileup\" \"path/to/bam/files/*.bam\" call-set -S (for SOLiD)"
	exit 0
end

pileup_files = Dir.glob(ARGV[0])
bam_files = Dir.glob(ARGV[1])
# #bam_files = ARGV[1]
STDERR.puts "Starting up..."
tn_gen = True_negative_site_generator.new(pileup_files, bam_files, ARGV[2], ARGV[3]=='-S')
STDERR.puts "running..."
tn_gen.run()
STDERR.puts "DONE"
