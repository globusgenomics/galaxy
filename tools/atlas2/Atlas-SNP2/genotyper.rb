#!/usr/bin/env ruby

usage='''
Program: v  generate vcf format output of a sample based on Atlas-SNP2 results file
Usage:   vcf.rb [Atlas-SNP2 result file] [sample name] [postrior probability cutoff] [minimal coverage] > [outputfile]

'''

if ARGV.size < 4
  puts usage
  exit
end

def passes_strandedness_test(reads_info, num_alt_reads,total_cov )
  return true if total_cov.to_i < 16 #with so few total reads, it is quite possible they are all one direction
  plus_dir = reads_info.scan(/\)([+-])[acgtACGTN]/).count(["+"])
  minus_dir = num_alt_reads.to_i - plus_dir
  return false if [plus_dir.to_f/num_alt_reads.to_f, minus_dir.to_f/num_alt_reads.to_f].min < 0.01
  return true
end



atlas_results = open(ARGV[0],'r')
sample_name = ARGV[1]
cutoff = ARGV[2].to_f
min_coverage = ARGV[3].to_f

puts '##fileformat=VCFv4.0'
puts "##fileDate=#{Time.now.strftime("%Y%m%d")}"
puts "##source=Atlas-Indel2-Illum-Exome"
puts "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
puts "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
puts "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
puts "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">"
puts "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major Variant Read Depth\">"
puts "##FILTER=<ID=low_VariantReads,Description=\"Number of variant reads is less than 3\">"
puts "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than 0.1\">"
puts "##FILTER=<ID=single_strand,Description=\"All variant reads are in a single strand direction\">"
puts "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than #{min_coverage.to_i}\">"
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t#{sample_name}\n"

atlas_results.each do |line|
  next if line =~ /^refName/
  cols = line.split("\t")
  ref, pos, ref_allele, var_allele, var_num, ref_num, pr, reads_info = cols[0].sub("chr",""), cols[1], cols[2], cols[3], cols[5].to_i, cols[7].to_i, cols[15].to_f, cols[18]
  n = ref_num + var_num
  
  next if pr < cutoff

  filter = "" 

  if n < min_coverage
    filter += "low_coverage;"
  end
  
  if not passes_strandedness_test(reads_info, var_num, n)  
    filter += "single_strand;"
  end

  if var_num < 3
    filter += "low_VariantReads;"
  end

  snp_qual = (-10 * Math.log10(1-pr+0.000001)).round
  
  if var_num.to_f/n <= 0.1
    genotype = "0/0"
    filter += "low_VariantRatio"
  elsif var_num.to_f/n > 0.1 and var_num.to_f/n < 0.9
    genotype = "0/1"
  else
    genotype = "1/1"
  end
  
  if filter == ""
    filter = "PASS" 
  else
    filter = filter.chomp(";")
  end

#  if genotype != "0/0"
    print ref,"\t",pos,"\t",".","\t",ref_allele,"\t",var_allele,"\t",snp_qual,"\t",filter,"\t",".","\t","GT:VR:RR:DP:GQ","\t","#{genotype}:#{var_num}:#{ref_num}:#{n}:.\n"
#  end

end
