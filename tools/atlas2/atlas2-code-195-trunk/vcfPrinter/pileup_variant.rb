 $:.unshift File.join(File.dirname(__FILE__),'.')

require 'pileupSiteInfo.rb'

class PileupLine

    attr_reader :chr, :coor, :refBase, :pos, :format, :total_depth

    def initialize(line)
        @values = line.strip.split("\t")
        @chr = @values[0]
        @coor = @values[1]
        @refBase = @values[2]
        @pos = "#{@chr}_#{@coor}"
        @format = "GT:VR:RR:DP:GQ:FT"
        @total_depth = @values[3].to_i
    end

    def ref_depth_snp
        siteObj=SiteInfo.new(@values[4])
        return siteObj.reference_depth_snp
    end

    def ref_depth_indel
        siteObj=SiteInfo.new(@values[4])
        return siteObj.reference_depth_indel
    end

    def genotype
        return './.'
    end

    def indels
        indels=@values[4].scan(/[+|-][0-9]+[ACGTNacgtn]+/)
        return indels.uniq
    end

    def var_depth
        return '.'
    end

    def genotype_qual
        return '.'
    end

    def genotype_filter
        return '.'
    end

    def genotypeInfo(varType)
        if varType=='SNP'
            return "#{self.genotype}:#{self.var_depth}:#{self.ref_depth_snp}:#{@total_depth}:#{self.genotype_qual}:#{self.genotype_filter}"
        else
            return "#{self.genotype}:#{self.var_depth}:#{self.ref_depth_indel}:#{@total_depth}:#{self.genotype_qual}:#{self.genotype_filter}"
        end
    end

end
