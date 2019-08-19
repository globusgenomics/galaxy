# usage ruby download_blah.rb sites file.gz 
header = false
File.open(ARGV[0], 'r').each_line do |line|
	next if line[0]=='#'
	cols=line.split("\t")
	chr = cols[0]
	coor = cols[1]
	if(header)
		h=""
	else
		h="-h"
		header= true
	end
	#puts `tabix #{h} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr#{chr}.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz #{chr}:#{coor}-#{coor}`.chomp
	puts `tabix #{h} ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110816_lc_indels_gls/ALL.chr#{chr}.VQSR_v0_gls_of_only_passing_sites.20101123.indels.genotypes.vcf.gz #{chr}:#{coor}-#{coor}`.chomp
	#puts `tabix #{h} #{ARGV[1]} #{chr}:#{coor}-#{coor}`.chomp
end
