$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 2)
	puts "USAGE: ruby intersect_vcfs.rb vcf_file1 vcf_file2 [0 (window size)]\ndoes not play nice with X and Y chromosomes"
	exit 0
end

def get_line(file)
	begin
		return file.readline.chomp
	rescue EOFError
		return :EOF
	end
end

def allele_matches(vcf_line1, vcf_line2)
	vcf_line1.alt.each do |alt1|
		alt1 = alt1.to_s
		vcf_line2.alt.each do |alt2|
			alt2 = alt2.to_s
			ref1 = vcf_line1.ref.to_s
			ref2 = vcf_line2.ref.to_s
			if(alt1.to_s.length > 1 && ref1.to_s.length > 1)
				ref1,alt1 = simplify(ref1.to_s,alt1.to_s)
			end
			if(alt2.to_s.length > 1 && ref2.to_s.length > 1)
				ref2,alt2 = simplify(ref2.to_s,alt2.to_s)
			end
			if(ref1 == ref2 && alt1 == alt2)
				return true
			end
		end
	end
	return false
end


def simplify(ref, alt)
	if(ref.length < alt.length)
		extra_length = ref.length - 1
		new_ref = ref[0]
		new_alt = alt[0..(alt.length-extra_length)]
	elsif(ref.length > alt.length)
		extra_length = alt.length - 1
		new_alt = alt[0]
		new_ref = ref[0..(ref.length-extra_length)]
	else
		raise "non-indel site!"
	end
	return [new_ref, new_alt]	
end


window = ARGV[1].to_i
labels1=nil
labels2=nil
window = 0
if(!ARGV[2].nil?)
	window = ARGV[2].to_i
end
prev_line = nil	
file1 = File.open(ARGV[0], 'r')
file2 = File.open(ARGV[1], 'r')
writer1=File.open("#{ARGV[0]}.inter#{window}", 'w')
writer2=File.open("#{ARGV[1]}.inter#{window}", 'w')
line1 = get_line(file1)
line2 = get_line(file2)
file2_buffer = Array.new
while(line1 != :EOF && (line2 != :EOF || file2_buffer.length > 0))
	if(line1[0,1] == '#')
		if(line1[0,6] == '#CHROM')
			labels1 = line1.split("\t")
		end
		writer1.puts(line1)
		line1 = get_line(file1)
		next
	end
	if(line2 != :EOF && line2[0,1] == '#')
		if(line2[0,6] == '#CHROM')
			labels2 = line2.split("\t")
			writer2.puts(line2)
			line2 = get_line(file2)
			vcf_line2 = Vcf_line.read_line(line2, false, labels2)
			next
		end
		writer2.puts(line2)
		line2 = get_line(file2)
		next
	end
	vcf_line1 = Vcf_line.read_line(line1, false, labels1)
	new_buffer = Array.new
	file2_buffer.each do |vcf_line| # remove sites no longer within range from the buffer
		if(vcf_line1.chr == vcf_line.chr && (vcf_line1.pos - vcf_line.pos).abs <= window)
			new_buffer.push vcf_line
		end
	end
	file2_buffer = new_buffer	
	while(line2 != :EOF && (vcf_line1.chr.to_i > vcf_line2.chr.to_i || (vcf_line1.chr == vcf_line2.chr && vcf_line1.pos -  vcf_line2.pos > window))) # skip file2 lines as needed
		line2 = get_line(file2)
		vcf_line2 = Vcf_line.read_line(line2, false, labels2) unless line2 == :EOF
	end
	while(line2 != :EOF && vcf_line1.chr == vcf_line2.chr && (vcf_line1.pos - vcf_line2.pos).abs <= window) # add new sites within range to the buffer
		file2_buffer.push(vcf_line2)
		line2 = get_line(file2)
		vcf_line2 = Vcf_line.read_line(line2, false, labels2) unless line2 == :EOF
	end

	file2_buffer.each do |vcf_line2|
		if(allele_matches(vcf_line1, vcf_line2))
			writer1.puts vcf_line1.print
			writer2.puts vcf_line2.print
			break
		end
	end
	line1 = get_line(file1)
	
end
