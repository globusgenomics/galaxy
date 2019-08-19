#!/usr/bin/ruby
$:.unshift File.join(File.dirname(__FILE__),'.')
require 'vcf_line'

if(ARGV.length < 2)
	puts "USAGE: ruby merge_duplicate_events.rb vcf_file window_size -s (merge_samples)"
	exit 0
end

window = ARGV[1].to_i
merge_samples = true if ARGV[2] == '-s'
labels=nil
prev_line = nil	
merge_count = 0
sample_replace_count = 0
warnings = 0
File.open(ARGV[0], 'r').each_line do |line|
	if(line[0,1] == '#')
		if(line[0,6] == '#CHROM')
			labels = line.split("\t")
		end
		puts line
		next
	end
	vcf_line = Vcf_line.read_line(line, false, labels)
	if(prev_line)
		chr = vcf_line.chr
		pos = vcf_line.pos
		if(chr == prev_line.chr && (prev_line.pos - pos).abs <= window && prev_line.ref == vcf_line.ref && (prev_line.alt & vcf_line.alt).length > 0) # is a match
			merge_count += 1
			vcf_line_to_keep = nil
			vcf_line_to_drop = nil
			ac=vcf_line.info[:"AC"]
			ac=ac.max if ac.class == Array
			prev_ac=prev_line.info[:"AC"]
			prev_ac=prev_ac.max if prev_ac.class == Array
			if(vcf_line.qual > prev_line.qual || (ac != nil && vcf_line.qual == prev_line.qual && ac >= prev_ac))
				vcf_line_to_keep = vcf_line
				vcf_line_to_drop = prev_line
			else
				vcf_line_to_keep = prev_line
				vcf_line_to_drop = vcf_line
			end
			if(merge_samples)
				vcf_line.sample_names.each do |name|
					keep_sample = vcf_line_to_keep.samples[name]
					drop_sample = vcf_line_to_drop.samples[name]
					# only look to merge sample if the vcf_line to keep is non-variant
					if(keep_sample[:"GT"].to_s =~ /[1-9]/)
						if(drop_sample[:"GT"].to_s =~ /[1-9]/)
							STDERR.puts "WARNING: merging lines #{chr}:#{prev_line.pos}-#{pos} which have nonvariant genotypes at both sites in sample #{name}"
							warnings += 1
						end
					else
						if(drop_sample[:"GT"].to_s =~ /[1-9]/)
							vcf_line_to_keep.samples[name]=drop_sample
							sample_replace_count += 1
						end
					end
				end
			end
			# clean up info that may not be accurate now
			vcf_line_to_keep.info.delete(:"AC")
			vcf_line_to_keep.info.delete(:"AN")
			vcf_line_to_keep.info.delete(:"NS")
			vcf_line_to_keep.info.delete(:"AAF")
			prev_line = vcf_line_to_keep
		else
			puts prev_line
			prev_line = vcf_line
		end
	else
		prev_line = vcf_line
	end
end
puts prev_line.print() unless prev_line.nil?
