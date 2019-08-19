
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'bed_file.rb'

if(ARGV.length < 2)
	puts "USAGE: ruby bed_filter.rb tab_delimited_file bed_file"
	exit 0
end

labels=nil
bed = Bed_file.new(ARGV[1])
	
File.open(ARGV[0], 'r').each_line do |line|
	if(line[0,1] == '#')
		puts line
		next
	end
	line.chomp!
	cols = line.split("\t")
	chr = cols[0]
	chr.slice!('chr')
	if(bed.pos_included?(chr.to_sym, cols[1].to_i))
		puts line
	end
end
		

