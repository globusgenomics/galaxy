
class Locus
	attr_reader :chr, :coor
	
	def initialize(chr, coor)
		chr=chr.to_s
		if(@chr=~/\d+/)
			@chr=chr.to_i
		else
			@chr=chr
		end
		@coor=coor.to_i
	end
	
	
	def <=>(other)
		case chr <=> other.chr
		when 0
			return coor <=> other.coor
		when -1
			return -1
		when 1
			return 1
		when nil
			if(chr.class == String && other.chr.class == Fixnum)
				return 1
			else
				return -1
			end
		else
			raise "invalid <=> result for locus chromosomes"
		end
	end

	def ==(other)
		return nil if other.nil?
		return chr == other.chr && coor == other.coor
	end


	def eql?(other)
		return self == other
	end

	def hash
		@chr.hash ^ @coor.hash
	end

	def to_s
		return "#{@chr}:#{@coor}"
	end

	def add(num)
		return Locus.new(@chr, @coor+num)
	end

	def subtract(num)
		return Locus.new(@chr, @coor-num)
	end
end
