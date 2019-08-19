
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'
require 'bed_file.rb'

if(ARGV.length < 2)
	puts "USAGE: ruby select_indels.rb vcf_file frameshift|inframe|1bp|insertion|deletion|1bpFilt|P|bed|AC|INFO|ACQual|MAF|AF|singSamples|siteList|length|notSiteList arg(if needed)"
	exit 0
end

select = ARGV[1]
arg = ARGV[2]
labels=nil
if(select == 'bed')
	bed = Bed_file.new(arg)
end
	
File.open(ARGV[0], 'r').each_line do |line|
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
		end
		puts line
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	#if(vcf_line.alt.length > 1)
	#	STDERR.puts "skipping multi-allelic site: #{line}"
	#	next
	#end
	0.upto(vcf_line.alt.length-1) do |i|
		case select
		when "frameshift"
			if((vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length) % 3 != 0 && (vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length) != 0)
				puts line
				break
			end
		when "inframe"
			if((vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length) % 3 == 0 && (vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length) != 0)
				puts line
				break
			end
		when "1bp"
			if( (vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length).abs == 1 )
				puts line
				break
			end
		when "insertion"
			if( vcf_line.alt[i].to_s.length > vcf_line.ref.to_s.length )
				puts line
				break
			end
		when "deletion"
			if( vcf_line.alt[i].to_s.length < vcf_line.ref.to_s.length )
				puts line
				break
			end
		when "1bpFilt"
			# if((vcf_line.alt[0].to_s.length - vcf_line.ref.to_s.length) != -1 || vcf_line.info[:P] >= arg.to_f) # filter 1bp deletions only
			if(((vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length) != -1 && (vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length) != 1)  || vcf_line.info[:P] >= arg.to_f) #filt 1bp insertions and deletions
				puts line
				break
			end
		when "P"
			if(vcf_line.info[:P] >= arg.to_f)
				puts line
				break
			end
		when 'bed'
			if(bed.pos_included?(vcf_line.chr.to_sym, vcf_line.pos))
				puts line
				break
			end
		when 'AC'
			low,high = arg.split(',')
			low = low.to_i
			high = high.to_i
			if(vcf_line.info[:"AC"].class == Array)
				sum = 0
				vcf_line.info[:"AC"].each {|x| sum += x.to_i}
				if(sum>low && sum<=high)
					puts line
					break
				end
			elsif(vcf_line.info[:"AC"]>low && vcf_line.info[:"AC"]<=high)
				puts line
				break
			end
		when 'INFO'
			id,low,high = arg.split(',')
			low = low.to_f
			high = high.to_f
			if(vcf_line.info[id.to_sym]>=low && vcf_line.info[id.to_sym]<=high)
				puts line
				break
			end
		when 'ACQual' #filter lines with an AC>= the specified count and with a QUAL value less than the specified cutoff
			args_cols = arg.split(',')
			ac = args_cols[1].to_i
			cutoff = args_cols[1].to_f
			if(vcf_line.info[:"AC"]>ac || vcf_line.qual >= cutoff)
				puts line
				break
			end
		when 'MAF'
			low,high = arg.split(',')
			low = low.to_f
			high = high.to_f
			maf = vcf_line.info[:"MAF"]
			if(maf.class == Array)
				sum = 0
				maf.each {|x| sum += x.to_f}
				maf = sum
			end
			if(maf>=low && maf<high)
				puts line
				break
			end
		when 'AF'
			low,high = arg.split(',')
			low = low.to_f
			high = high.to_f
			maf = vcf_line.info[:"AF"]
			if(maf.class == Array)
				sum = 0
				maf.each {|x| sum += x.to_f}
				maf = sum
			end
			if(maf>=low && maf<high)
				puts line
				break
			end
		when 'AAF'
			low,high = arg.split(',')
			low = low.to_f
			high = high.to_f
			maf = vcf_line.info[:"AAF"]
			if(maf.class == Array)
				sum = 0
				maf.each {|x| sum += x.to_f}
				maf = sum
			end
			if(maf>=low && maf<high)
				puts line
				break
			end
		when 'singSamples'
			if(vcf_line.info[:"AC"]==1)
				vcf_line.samples.each_pair do |name, sample|
					if(sample[:"GT"] =~/1/)
						puts "#{vcf_line.chr}:#{vcf_line.pos}	#{name}	#{vcf_line.ref}	#{vcf_line.alt.join(',')}"
						break
					end
				end
			end
		when 'siteList'
			if(`grep -cP "^#{vcf_line.chr}	#{vcf_line.pos}(	|$)" #{arg}`.chomp.to_i > 0)
				puts line
				break
			end
		when 'notSiteList'
			if(`grep -cP "^#{vcf_line.chr}	#{vcf_line.pos}(	|$)" #{arg}`.chomp.to_i == 0)
				puts line
				break
			end
		when "length"
			low,high = arg.split(',')
			low = low.to_f
			high = high.to_f
			length = (vcf_line.alt[i].to_s.length - vcf_line.ref.to_s.length).abs
			if(length>=low && length<high)
				puts line
				break
			end
		else
			STDERR.puts "I don't know how to select indels by #{select}, quitting"
			exit 1
		end
	end
end
		

