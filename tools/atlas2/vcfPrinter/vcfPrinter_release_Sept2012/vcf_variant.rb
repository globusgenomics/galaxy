
class VcfLine
	attr_reader :chr, :coor, :pos, :format, :uniqPos, :sampleInfo, :qual
	attr_accessor :refBase, :altBase, :filter

	def initialize(line)
		@values = line.strip.split("\t")
		@chr = @values[0]
		@coor = @values[1]
		@refBase = @values[3].upcase
		@altBase = @values[4].upcase 
		@uniqPos = "#{@chr}_#{@coor}_#{@refBase}_#{@altBase}"
		@pos = "#{@chr}_#{@coor}"
		@format = @values[8]   #"GT:VR:RR:DP:GQ:FT"
		@qual = @values[5].to_i
	end

	def updateQual(qualScore)
        #while vcfPrinter is looping through variants from various samples, qual score is updated if this method finds a qual score greater than present score.
        if qualScore == '.'

        else
                if @qual >= qualScore.to_i

                else
                        @qual=qualScore
                end
        end
        end
	

	def filter
		return @values[6]
	end

	def info
		return @values[7]
	end

	def id
		return @values[2]
	end

	def genotype
	#If the filter field is anything other than PASS then the genotype is changed to ./.
		if self.filter == 'PASS'
			return @values[-1].split(":")[0]
		else
			return './.'
		end
	end

	def var_depth
		return @values[-1].split(":")[1]
	end

	def ref_depth
		return @values[-1].split(":")[2]
	end

	def total_depth
		return @values[-1].split(":")[3]
	end

	def genotype_qual
		return @values[-1].split(":")[4]
	end

	def genotype_filter
		return self.filter
	end

	def type
                if @values[3].length > @values[4].length
                        return 'Del'
                elsif @values[3].length < @values[4].length
                        return 'Ins'
                else
                        return 'SNP'
                end
        end

	def genotypeInfo
		return "#{self.genotype}:#{self.var_depth}:#{self.ref_depth}:#{self.total_depth}:#{self.genotype_qual}:#{self.genotype_filter}"
	end

	def sampleInfoCols
		return @values[9..-1]
	end

	def updateSampleInfo(genoInformation)
	#Links the alt alleles with samples they came from
		#define an new array for sample genotype information
		@sampleInfo=Array.new()
		#get the allele number by counting the number of allele listed in altBase column
		alleleNum=@altBase.strip().split(',').length
		#get to be replaced sample geno info from vcfline
		oldSampleInfo=sampleInfoCols()
		#loop through them
                oldSampleInfo.each_with_index do |varInfo, index|
			#if the sample geno info is found in noChange bin then dont change anything and add it to the new sample geno list
			if !genoInformation[:noChange][index].nil?
				@sampleInfo[index]=genoInformation[:noChange][index]
			#if it is found in the change bin then change the genotype from 1/1 to new allele number 4/4
			elsif !genoInformation[:change][varInfo].nil?
				@sampleInfo[index]="#{varInfo.split(':')[0].gsub('1',alleleNum.to_s)}:#{varInfo.split(':')[1..-1].join(':')}"
			#all the other add to new list
			else
				@sampleInfo[index]=varInfo
			end
		end
        end

	def update(vcfLineObj)
	#This function is used when collapsing variants. variables to be updated refBase,altBase,Filter
		vcfLineObj.altBase.split(',').each do |vcfAltAllele|
			#if @altBase.split(',').include?(vcfAltAllele)
			#else
				@altBase="#{vcfAltAllele},#{@altBase}"
			#end
		end
		vcfLineObj.refBase.split(',').each do |vcfRefAllele|
			#if @refBase.split(',').include?(vcfRefAllele)
                	#else
                        	@refBase="#{vcfRefAllele},#{@refBase}"
	                #end
		end
		if self.filter.split(';').include?(vcfLineObj.filter)
                else
                        @filter="#{vcfLineObj.filter},#{@filter}"
                end
	end

end

#Testing vcf.rb
#v = VcfLine.new("1	113054374	.	CTTG	C	23	PASS	AC=2;AN=4;DP=4600;NS=65	GT:VR:RR:DP:GQ	./.:.:.:128:.")
#puts v.genotypeInfo

