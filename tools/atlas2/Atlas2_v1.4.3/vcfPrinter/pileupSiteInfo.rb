#!/usr/bin/python

class SiteInfo

    def initialize(pileupSiteInfoLine)
        @siteInfoLine=pileupSiteInfoLine
    end

    def reference_depth_snp
        return @siteInfoLine.count(".,").to_i
    end 

    def reference_depth_indel
        vars=@siteInfoLine.scan(/([\w|\W])([+|-])([0-9])+[ACGTNacgtn]/)
        totalRef = @siteInfoLine.count(".,").to_i
        return totalRef if vars.length == 0
        vars.each do |var|
            totalRef=totalRef - 1 if [',', '.'].include?(var[0])
        end
        return totalRef
    end

    def get_indel_alleles
        vars=@siteInfoLine.scan(/([\w|\W][+|-][0-9]+[ACGTNacgtn]+)/)
        return vars.uniq
    end

    def get_indel_allele_counts
        alleleCount=Hash.new
        vars=@siteInfoLine.scan(/([\w|\W][+|-][0-9]+[ACGTNacgtn]+)/)
        vars.each do |var|
            if alleleCount.include?(var)
              alleleCount[var]+=1
            else
              alleleCount[var]=1
            end
        end
        return alleleCount
    end

    def get_major_indel_allele
    #In Beta
        indels=self.get_indel_allele_counts
        majorAllele=String.new()
        max=0
        indels.each do |allele, count|
            if count > max 
              majorAllele = allele
              max = count
            else
              next
            end
        end
        return majorAllele
    end
end 
