#!/usr/bin/ruby
if(ARGV.length < 4)
	puts "USAGE: ruby parellele_process_big_vcf.rb [file] [name_prefix] [num_jobs] [ram] \"command %%% > %%%.result]\""
	exit 0
end
file = ARGV[0]
name=ARGV[1]
num_jobs = ARGV[2].to_i
print "Sizing up the file..."
lines_per_job=`wc -l #{file}`.chomp.to_i / num_jobs
puts "DONE"
ram=ARGV[3]
# always enter file as $file.vcf
# command="ruby /stornext/snfs6/1000GENOMES/challis/lib/VCF-Utils/add_allele_freq.rb $file.vcf -clean | grep -v '^#'> $file.vcf.clean"
command=ARGV[4].dup
command.gsub!(/%%%/,'$file.vcf')
print "Splitting file..."
`split -l #{lines_per_job} #{file} #{name}.#{$$}.`
puts "DONE"
`head -n 500 #{file} | grep "^#" > tmp.header.#{$$}`
print "Submitting jobs..."
`for file in #{name}.#{$$}.* ; do echo "cat tmp.header.#{$$} > $file.vcf; grep -v '#' $file >> $file.vcf; #{command}" | msub -q analysis -d $PWD -l pmem=#{ram} -N #{name}.$file; done`
puts "DONE"
puts "Mapping done, you have to reduce it yourself. Remember to strip out the headers!"

