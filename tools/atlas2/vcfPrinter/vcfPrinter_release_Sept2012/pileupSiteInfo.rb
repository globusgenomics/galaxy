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

end 
