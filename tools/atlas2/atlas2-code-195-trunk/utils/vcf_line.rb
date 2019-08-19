class Vcf_line
	include Comparable
	require 'Locus'
	attr_accessor  :locus, :id, :ref, :alt, :qual, :filter, :info, :format, :samples, :source_file, :sample_names


	def initialize(chr, pos, ref, alt, format, samples=Hash.new)
		@locus = Locus.new(chr,pos)
		@ref = ref.to_sym
		@alt = alt.split(',')
		@alt.collect! {|base| base.to_sym}
		if(format.nil?)
			@format = Array.new
		else
			@format = format.split(":")
			@format.collect! {|value| value.to_sym}
		end
		@info = Hash.new
		@samples=samples
	end

 
       # parse a VCF line, 
	# line_str: the vcf line as a string
	# filter: bool, should lines with a filter raise exceptions?
	# labels: the head labels line as an array of strings [line.split("\t")]	
	# source_file: optional, the path (or is it name?) to the source file
	def self.read_line(line_str, filter, labels, source_file=nil)
		begin
			line_str.gsub!(/\s+$/, "")
			line = line_str.split("\t")
			if(filter && line[6] != 'PASS' && line[6] != '.')
				raise "filter"
			end
			vcf_line = Vcf_line.new(line[0], line[1], line[3], line[4], line[8])
			vcf_line.id = line[2].to_sym
			vcf_line.qual = line[5].to_i
			vcf_line.filter = line[6].to_sym
			vcf_line.parse_info(line[7])
			vcf_line.sample_names = Array.new
			if line.length != labels.length
				#STDERR.puts "WARNING: line #{vcf_line.chr}:#{vcf_line.pos} size #{line.length} does not match labels size #{labels.length}" 
			else
				(9...line.size).each do |index|
	     				sample_name = labels[index].chomp
					vcf_line.sample_names.push(labels[index].chomp)
					vcf_line.parse_sample(line[index], sample_name)
				end
			end
			#    rescue => detail
			#      puts "ERROR: One or more lines in the vcf file is in the incorrect format"
			#      puts detail.message
			#      Process.exit 0
			vcf_line.source_file=source_file # a reference back to the Vcf_file object this line belongs to (if any)
			return vcf_line
		end
	end #initialize
	
	def chr
		return @locus.chr
	end
	
	def pos
		return @locus.coor
	end


	def print
		return "#{chr}\t#{pos}\t#{@id}\t#{@ref}\t#{@alt.join(',')}\t#{@qual}\t#{@filter}\t#{print_info}\t#{@format.join(':')}\t#{print_samples}".chomp.chomp
	end

	def to_s
		print
	end

	def print_bin(conf)
		return "#{chr}\t#{pos}\t#{@ref}#{@alt[0].to_s.downcase}\t#{print_samples_bin(conf)}"
	end

	def add_filter(filter)
		if(@filter == :"PASS" || @filter == :".")
			@filter = filter.to_sym
		else
			@filter = "#{@filter};#{filter}".to_sym
		end
	end

	def is_frameshift(alt_i=0)
		return indel_length(alt_i) % 3 != 0	
	end

	def indel_length(alt_i=0)
		return @alt[alt_i].to_s.length - @ref.to_s.length
	end

	def get_pileup_alt(alt_i=0)
		length = indel_length(alt_i)
		ref, alt = Vcf_line.simplify(@ref.to_s, @alt[alt_i].to_s)
		if(ref[0] == alt[0])
			anchor=''
		else # complex
			anchor = alt[0]
		end
		if(length > 0) # insertion 
			alt[0]=''
			return "#{anchor}+#{length}#{alt}"
		elsif(length < 0) # deletion
			ref[0]=''
			return "#{anchor}#{length}#{ref}"
		else
			raise "this method is only designed for indels"
		end
	end

	def <=> (partner)
		return @locus <=> partner.locus
	end #<=>

	def == (partner)
		if( (self <=> partner) == 0 )
			if(@ref == partner.ref && @alt == partner.alt)
				return true
			end
		end
		return false
	end


	def get_chrom_num
		if(chr.to_s =~ /^(\d+)(.*)$/)
			if($2 == nil || $2 == '')
				return $1.to_i.to_f
			end
			return ($1.to_i.to_f ) + (1.0 / $2.hash.to_f.abs )
		end
		return 1000.0 + chr.to_s.hash.to_f.abs  #this will give this chromosome a unique id above 10002
	end #get_chrom_num


	#return true if the indels are roughly equivalent and should be merged
	def is_equivalent(other, fuzz=5)
		if(chr == other.chr)
			if( (pos - other.pos).abs < (fuzz + 1) ) # coor must be no more than fuzz bp apart
				#indel length must be no more than 2 bp different
				if(@ref.to_s.length==1 && other.ref.to_s.length==1 && (@alt[0].to_s.length - other.alt[0].to_s.length).abs < 3)
					return true
				elsif(@alt[0].to_s.length==1 && other.alt[0].to_s.length==1 && (@ref.to_s.length - other.ref.to_s.length).abs < 3)
					return true
				end
			end
		end
		return false
	end
	
	
	def self.simplify(ref, alt)
		if(ref.length < alt.length)
			extra_length = ref.length
			new_ref = ref[0]
			new_alt = alt[0..(alt.length-extra_length)]
		elsif(ref.length > alt.length)
			extra_length = alt.length
			new_alt = alt[0]
			new_ref = ref[0..(ref.length-extra_length)]
		else
			return [ref, alt]
		end
		return [new_ref, new_alt]
	end


	def parse_info(info_str)
		@info = Hash.new
		return if(info_str == nil)
		info_arr = info_str.split(';')
		info_arr.each do |entry_str|
			entry = entry_str.split('=')
			key = entry[0].to_sym
			if(entry[1] =~ /,/)
				values = entry[1].split(',')
				values.each {|value| parse_value(value)}
				@info[key] = values
			elsif(entry[1] == nil)
				@info[key] = nil
			else
				value = parse_value(entry[1])
				@info[key] = value
			end
		end
	end

	def parse_sample(sample_str, sample_name)
		sample = Hash.new
		if(sample_str == './.' || sample_str == '.')
			@samples[sample_name]=sample
			return
		end
		sample_arr = sample_str.split(':')
		if(sample_arr.size != @format.size)
			raise "#{chr}:#{pos} The sample #{sample_name} format #{sample_arr.size} does not match the specified format #{@format.size} -- #{sample_arr}"
		end
		end_index = sample_arr.size-1
		(0..end_index).each do |index|
			value = sample_arr[index]
			parse_value(value)
			sample[@format[index]] = value
		end
		@samples[sample_name]=sample
	end



	def print_sample(sample_name, str='')
		sample = @samples[sample_name]
		# these lines were added for a specialized imputation variant and should not be uncommented
		# str = "#{str}#{sample[:"RR"]} #{sample[:"VR"]}\t"
		# next
		@format.each do |key|
			if(sample.nil? || sample.length == 0)
				#STDERR.puts "Warning: sample #{sample_name} was not found in the sample hash" if sample.nil?

				str = "#{str}.:"
			elsif(sample[key].instance_of?(Array))
				str = "#{str}#{sample[key].join(',')}:"
			elsif(!sample.include?(key))
				str = "#{str}.:"
			else
				str = "#{str}#{sample[key]}:"
			end
		end
		str.chop!
		return str
	end


	private

	

	def print_info
		if(@info.length == 0 || (@info.length == 1 && @info.include?(:".")))
			return '.'
		end
		@info.delete(:".")
		str = ''
		@info.each do |key, value|
			if(value == nil)
				str = "#{str}#{key};"
			else
				if(value.instance_of?(Array))
					str = "#{str}#{key}=#{value.join(',')};"
				else
					str = "#{str}#{key}=#{value};"
				end
			end
		end
		str.chop!
		return str
	end

	

	def print_samples
		str = ''
		if( !@sample_names.nil? && @sample_names.length > 0 && @sample_names[0] != nil )
			@sample_names.each do |sample_name|
				str=print_sample(sample_name, str)
				str = "#{str}\t"
			end
		end
		str.chop!
		return str
	end #print_samples





	def print_samples_bin(conf)
		str = ''
		if( @sample_names.length > 0 && @sample_names[0] != nil )
			@sample_names.each do |sample_name|
				sample = @samples[sample_name]
				if(sample == nil || sample.length==0 || sample[:"GT"] =~/\.[\\\/|]\./)
					str="#{str}0.34 0.33"
				elsif(!sample[:"GQ"].nil?) # use the genotype quality if it is there
					str="#{str}#{sample[:"GQ"].split(',')[0]} #{sample[:"GQ"].split(',')[1]}"
				else
					raise "missing sample GQ values!"
				end
				#elsif( sample[:"GT"] =~ /[123456789][\/\\|][123456789]/ )
				#	str="#{str}#{((100000.0*((1.0-conf)/2.0)).round.to_f)/100000.0} #{conf}"
				#elsif( sample[:"GT"] =~ /[123456789]/ ) # heterozygous
				#	str="#{str}#{conf} #{((100000*((1.0-conf)/2.0)).round.to_f)/100000.0}"
				#elsif( sample[:"GT"] =~ /0[\/\\|]0/ ) # homozygous
				#	str="#{str}#{((100000*((1.0-conf)/2.0)).round.to_f)/100000.0} #{((100000*((1.0-conf)/2.0)).round.to_f)/100000.0}"
				#else
				#	raise "unsupported case: #{sample}"
				#end
				str = "#{str}\t"
			end
		end
		str.chop!
		return str
	end #print_samples


	


	def parse_value(value)
		if(value =~ /^-*\d+$/)
			value = value.to_i
		elsif(value =~ /^-*\d+\.\d+$/)
			value = value.to_f
		else
			value = value.to_sym
		end
	end
	
	def equal_site(vcf_line)
		if(chr == vcf_line.chr && pos == vcf_line.pos)
			return true
		else
			return false
		end
	end

end #class


#test
#line_str = "20\t14370\trs6054257\tG\tA,T\t29\t0\tNS=3;DP=14;AF=0.5,0.2;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:-1,-1"
#line = Vcf_line.new(line_str)
#puts line_str
#puts line.print

