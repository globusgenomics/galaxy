$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_file.rb'

module Vcf_template

	def self.doit(filename)
		vcf_file = Vcf_file.new(filename)
		puts vcf_file.print_header()
		vcf_file.each_line do |line|
			#do something
			puts line.print
		end
	end

end



if __FILE__ == $0
	if(ARGV.length < 1 || ARGV[0] == "-h" || ARGV[0] == "--help")
		puts "USAGE: Template.rb vcf_file"
		exit 1
	end
	vcf_template.doit(ARGV[0])

end
