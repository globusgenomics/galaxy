class BedError < Exception
end
#require 'lib/serializer'
require 'read'


class Bed_file
	attr_reader :filename
	$bed_file_instance = nil

	def initialize(filename)
		@filename = filename
		@bed = Hash.new
    		File.open(filename, 'r').each_line do |line|
			cols = line.split("\t")
			cols[0].slice!('chr')
			chr = cols[0].to_sym
			@bed[chr] = [Array.new, Array.new] if @bed[chr] == nil
			@bed[chr][0].push( cols[1].to_i )
			@bed[chr][1].push( cols[2].to_i )
		end
		$bed_file_instance = self
	end

#	def self.load(filename)
#		$bed_file_instance = Serializer.load(filename)
#		raise "Error: passed file is not a serialized bed file" unless $bed_file_instance.class == Bed_file
#	end


	def self.instance
		raise "Bed_file.new must be called before Bed_file.instance" if $bed_file_instance == nil
		return $bed_file_instance
	end

	# returns true if the @bed has an entry that includes this position
	def pos_included?(chr, coor)
		return false if @bed[chr].nil?
		raise "Error: Invalid Bed File. It has a different number of start positions #{@bed[chr][0].length} than end positions #{@bed[chr][1].length}" if @bed[chr][0].length != @bed[chr][1].length
		return bin_search(coor, @bed[chr])
	end


	# returns true if the read is included, or nearly included in the bed file. IT WILL RETURN SOME READS THAT ARE NOT QUITE IN THE BED FILE, BUT NEARLY SO
	# it uses twice the read length as the read ending position in case of long deletions
#	def read_included?(read)
#		if(pos_included?(read.ref, read.pos) || pos_included?(read.ref, read.pos+(read.seq.length)))
#			return true
#		end
#		return false
#	end




	private

	# binary search helper for pos_included?
	def bin_search(search_pos, arr, start=0, stop=(arr[0].length-1))
		if stop - start < 3
			(start..stop).each do |i|
				if(search_pos > arr[0][i] && search_pos <= arr[1][i])
					return true
				end
			end
			return false
		end
		x = start + ((stop-start)/2)
		c = search_pos <=> arr[0][x]
		case c
		when 0
#			return true
			return bin_search(search_pos, arr, start, x-1)
		when -1
			return bin_search(search_pos, arr, start, x-1)
		else # 1
			return bin_search(search_pos, arr, x, stop)
		end
	end
end

# #MemoryProfiler.start #bed = Bed_file.new(ARGV[0]) #Bed_file.dump(bed)

