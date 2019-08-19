
class VcfLine
    attr_reader :chr, :coor, :pos, :uniqPos, :sampleInfo
    attr_accessor :refBase, :altBase, :qual, :filter

    def initialize(line)
        @values = line.strip.split("\t")
        @chr = @values[0].gsub('chr','')
        @coor = @values[1]
        @refBase = @values[3]
        @altBase = @values[4]
        @uniqPos = "#{@chr}_#{@coor}_#{@refBase}"
        @pos = "#{@chr}_#{@coor}"
        @qual = @values[5].to_i
        @filter = @values[6]
    end

    def updateQual(qualScore)
      #method to update qual score when looping thru samples
      if qualScore == '.'
        
      else

        if @qual >= qualScore.to_i

        else

          @qual=qualScore.to_i
        end
      end
    end

    def filter=(incoming_filter)
        if incoming_filter == '.' or @filter == 'PASS'

        else

            if incoming_filter == 'PASS'
              @filter = 'PASS'
            else
              incoming_filter.split(';').each do |x|
                @filter = @filter + ';' + x if !@filter.include? x
              end
            end
        end
    end


    def length
        #Function will return an array with length of all alt alleles
        len=Array.new
        alts=@altBase.strip.split(',')
        alts.each do |alt|
          len.push((@refBase.length-alt.length).abs)
        end
        return len
    end

    def format
      if @values[8].include? 'FT:AA'
        return "#{@values[8].strip}"
      else
        return "#{@values[8].strip}:FT:AA"
      end
    end

    def info
      return @values[7]
    end

    def id
      return @values[2]
    end

    def type
	if @values[3].length ==1 and @values[4].length ==1
	    return 'SNP'
        elsif @values[3].length > @values[4].length
            return 'Del'
        elsif @values[3].length < @values[4].length
            return 'Ins'
        else
            return 'ComplexVar'
        end
    end

    def genotypeInfo
      return "#{@values[9].strip}:#{self.filter}:#{self.altBase}"          
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
		@filter="#{vcfLineObj.filter};#{@filter}"
	end
    end

end

#Testing vcf.rb
#v = VcfLine.new("1    113054374    .    CTTG    C    23    PASS    AC=2;AN=4;DP=4600;NS=65    GT:VR:RR:DP:GQ    ./.:.:.:128:.")
#puts v.genotypeInfo

