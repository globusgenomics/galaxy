$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_editor.rb'

class Info_fixer < Vcf_editor

	def initialize(info_fields)
		super()
		@info_fields = info_fields
		@info_fields.map!{|key| key.to_sym}
		@allele_i = 0
	end

	def edit_line(line)
		if(line.locus == @prev_locus)
			@allele_i += 1
		else
			@allele_i = 0
		end
		@info_fields.each do |key|
			value = [line.info[key]]
			value.flatten!
			line.info[key]=value[@allele_i]
		end
		@prev_locus = line.locus
	end

end


if __FILE__ == $0
	if(ARGV.length < 1 ||ARGV[0] == "-h" || ARGV[0] == "--help")
		puts "USAGE: cat vcf_file | ruby fix_info_after_split.rb info,fields"
		exit 1
	end
	
	editor = Info_fixer.new(ARGV[0].split(","))
	editor.run
end
