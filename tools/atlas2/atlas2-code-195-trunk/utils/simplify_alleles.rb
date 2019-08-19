$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_file.rb'

module Simplify_alleles

	def self.doit()
		vcf_file = Vcf_file.new()
		puts vcf_file.print_header()
		vcf_file.each_line do |line|
			raise "This program only works for biallelic VCF lines" if line.alt.length > 1
			ref, alt = Vcf_line.simplify(line.ref.to_s, line.alt[0].to_s)
			line.ref = ref
			line.alt=[alt]
			puts line.print
		end
	end

end



if __FILE__ == $0
	if(ARGV[0] == "-h" || ARGV[0] == "--help")
		puts "USAGE: cat vcf_file | simplify_alleles.rb"
		exit 1
	end
	Simplify_alleles.doit()

end
