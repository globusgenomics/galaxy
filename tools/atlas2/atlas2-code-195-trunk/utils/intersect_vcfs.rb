$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 2)
	puts "USAGE: ruby intersect_vcfs.rb vcf_file1 vcf_file2\ndoes not play nice with X and Y chromosomes"
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
			begin
			if(alt1.to_s.length > 1 && ref1.to_s.length > 1)
				fref1,falt1 = simplify(ref1.to_s,alt1.to_s)
			else
				fref1 = ref1
				falt1 = alt1
			end
			if(alt2.to_s.length > 1 && ref2.to_s.length > 1)
				fref2,falt2 = simplify(ref2.to_s,alt2.to_s)
			else
				fref2 = ref2
				falt2 = alt2
			end
			rescue Exception => e
				puts e.message
				puts vcf_line1.print
				puts "###################################################################\n######################################"
				puts vcf_line2.print
				raise 'bad!'
			end
			if(fref1 == fref2 && falt1 == falt2)
				return true
			end
		end
	end
	return false
end


def simplify(ref, alt)
	if(ref.length < alt.length)
		extra_length = ref.length
		new_ref = ref[0]
		new_alt = alt[0..(alt.length-extra_length)]
	elsif(ref.length > alt.length)
		extra_length = alt.length
		new_alt = alt[0]
		new_ref = ref[0..(ref.length-extra_length)]
	else
		raise "non-indel site! #{ref}->#{alt}"
	end
	return [new_ref, new_alt]	
end


window = ARGV[1].to_i
labels1=nil
labels2=nil
prev_line = nil	
file1 = File.open(ARGV[0], 'r')
file2 = File.open(ARGV[1], 'r')
writer1=File.open("#{ARGV[0]}.inter", 'w')
writer2=File.open("#{ARGV[1]}.inter", 'w')
line1 = get_line(file1)
line2 = get_line(file2)
while(line1 != :EOF && line2 != :EOF)
	if(line1[0,1] == '#')
		if(line1[0,6] == '#CHROM')
			labels1 = line1.split("\t")
		end
		writer1.puts(line1)
		line1 = get_line(file1)
		next
	end
	if(line2[0,1] == '#')
		if(line2[0,6] == '#CHROM')
			labels2 = line2.split("\t")
		end
		writer2.puts(line2)
		line2 = get_line(file2)
		next
	end
	vcf_line1 = Vcf_line.read_line(line1, false, labels1)
	vcf_line2 = Vcf_line.read_line(line2, false, labels2)
	if(vcf_line1.chr == vcf_line2.chr && vcf_line1.pos == vcf_line2.pos)
		if(allele_matches(vcf_line1, vcf_line2))
			writer1.puts line1
			writer2.puts line2
		end
		line1 = get_line(file1)
		line2 = get_line(file2)
		next
	elsif(vcf_line1.chr.to_i > vcf_line2.chr.to_i || (vcf_line1.chr == vcf_line2.chr && vcf_line1.pos > vcf_line2.pos))
		line2 = get_line(file2)
		next
	elsif(vcf_line1.chr.to_i < vcf_line2.chr.to_i || (vcf_line1.chr == vcf_line2.chr && vcf_line1.pos < vcf_line2.pos))
		line1 = get_line(file1)
		next
	else
		raise "unhandled case!"
	end
	
end
