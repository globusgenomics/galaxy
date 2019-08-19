#!/hgsc_software/ruby/latest/bin/ruby
$:.unshift File.join(File.dirname(__FILE__),'.')
$: << Dir.pwd

raise RuntimeError, "settings.rb does not exist!" if !File.exist? ("#{Dir.pwd}/settings.rb")

require "#{Dir.pwd}/settings.rb"
require 'logger'
require 'getoptlong'
require 'tempfile'
require 'tmpdir'
require 'fileutils'
require 'vcf.rb'
require 'vcf_variant.rb'
require 'genotype.rb'
require 'clusterJobStats.rb'
require 'functions.rb'

Log = Logger.new("#{File.dirname(__FILE__)}/printer.log")
#Levels (DEBUG<INFO<WARN<ERROR<FATAL)
Log.level = Logger::INFO

class VcfOnCluster

	def initialize(cline, args)
          @njobs = args['jobs']
          #The chromosome size is based on HG19
          @chrs = Hash['1' => [249250621, @njobs], '2' => [243199373, @njobs], \
            '3' => [198022430, @njobs], '4' => [191154276, @njobs], \
            '5' => [180915260, @njobs], '6' => [171115067, @njobs], \
            '7' => [159138663, @njobs], '8' => [146364022, @njobs], \
            '9' => [141213431, @njobs], '10' => [135534747, @njobs],\
            '11' => [135006516, @njobs], '12' => [133851895, @njobs], \
            '13' => [115169878, @njobs], '14' => [107349540, @njobs], \
            '15' => [102531392, @njobs], '16' => [90354753, @njobs], \
            '17' => [81195210, @njobs], '18' => [78077248, @njobs], \
            '19' => [59128983, @njobs], '20' => [63025520, @njobs], 
            '21' => [48129895, @njobs], '22' => [51304566, @njobs], \
            'X' => [155270560, @njobs], 'Y' => [59373566, @njobs], 'MT' => [16569, @njobs]]
	  @outputFile=args['out_file']
          @single_chr = args['chr']
          @pass = args['pass']
          @hasPileup = args['hasPileup']
          @command_line = cline
	end

	def submitClusterJobs()
          chr_jobid = Array.new
	  chr_file = Array.new
          if @single_chr
            start_pos = 1
            chr = @single_chr.to_s
            chunk_size = @chrs[@single_chr][0] / @chrs[@single_chr][1]
            i=1
            until i > @chrs[@single_chr][1] do
              end_pos = i*chunk_size
              timeObj = Time.now
              tFileName = File.join(SCRATCH, "/#{timeObj.to_i.to_s}vcfjob_#{rand(1000000)}.vcf")
              randInt = rand(1000000)

              id=`echo \'#{@command_line} -o #{tFileName} -s #{chr}:#{start_pos}-#{end_pos}\' \
| qsub -q #{JOB_QUEUE} -l nodes=1:ppn=1,mem=6000mb \
-e #{SCRATCH}/#{timeObj.to_i.to_s}_vcfjob_chr_#{chr}_#{randInt}.err \
-o #{SCRATCH}/#{timeObj.to_i.to_s}_vcfjob_chr_#{chr}_#{randInt}.out -d $PWD`.gsub('.sug-moab', '')

              STDERR.puts "Submitted Chromosome #{chr}:#{start_pos}-#{end_pos} job to cluster jobid #{id.strip}"
              Log.info "Submitted Chromosome #{chr}:#{start_pos}-#{end_pos} job to cluster jobid #{id.strip}"

              chr_file.push(tFileName)
              chr_jobid.push(id.strip)
              start_pos = end_pos
              i+=1
            end
          else
            @chrs.each do |chr, value|
              start_pos = 1
              chunk_size = value[0] / value[1]
              i=1
              until i > value[1] do
                end_pos = i*chunk_size
                timeObj = Time.now
		tFileName = File.join(SCRATCH, "/#{timeObj.to_i.to_s}vcfjob_#{rand(1000000)}.vcf") 
                randInt = rand(1000000)

                id=`echo \'#{@command_line} -o #{tFileName} -s #{chr}:#{start_pos}-#{end_pos}\' \
| qsub -q #{JOB_QUEUE} -l nodes=1:ppn=1,mem=6000mb \
-e #{SCRATCH}/#{timeObj.to_i.to_s}_vcfjob_chr_#{chr}_#{randInt}.err \
-o #{SCRATCH}/#{timeObj.to_i.to_s}_vcfjob_chr_#{chr}_#{randInt}.out -d $PWD`.gsub('.sug-moab', '') 

                STDERR.puts "Submitted Chromosome #{chr}:#{start_pos}-#{end_pos} job to cluster jobid #{id.strip}"
                Log.info "Submitted Chromosome #{chr}:#{start_pos}-#{end_pos} job to cluster jobid #{id.strip}"

		chr_file.push(tFileName)
		chr_jobid.push(id.strip)
                start_pos = end_pos
                i+=1
              end
            end
	  end
	return chr_jobid, chr_file
	end

	def processJobs(chr_jobid, chr_file)
		while(trackJobs(chr_jobid))
			STDERR.puts "Processing"
		end
		chr_file.each_index do |fileIndex|
                        if !File.exists?(chr_file[fileIndex])
                          STDERR.puts "WARNING: Did not find #{chr_file[fileIndex]}"
                          Log.error "Did not find #{chr_file[fileIndex]}"

                          next
                        end
			if fileIndex==0
                          `cat #{chr_file[fileIndex]} >#{@outputFile}`
			else
			  `cat #{chr_file[fileIndex]} | grep -v '^#' >>#{@outputFile} `
			end
			#FileUtils.rm(outfile)
		end
		STDERR.puts "Done processing!"
	end

	def trackJobs(chr_jobid)
		while(jobsAreRunning(chr_jobid))
                  #STDERR.puts "Tracking jobs on cluster"
		  sleep(30)
		  return true
		end
	return false
	end

	def jobsAreRunning(chr_jobid)
		clusterJobTrackingObj = ClusterJobStatus.new
		jobsOnCluster = clusterJobTrackingObj.jobs
		#puts "Number of Jobs running on cluster: #{jobsOnCluster.length}"
		chr_jobid.each do |id|
			return true if jobsOnCluster.include?(id)	
		end
		return false
	end

	private :trackJobs, :jobsAreRunning

end

def singleCorePrinterJob(cline, out_file, onlyPass=false)
  if onlyPass
    `#{cline} -o #{out_file} -p`
  else
    `#{cline} -o #{out_file}`
  end
end

def quickPrinter(vcfs, out_file, onlyPass=false)
  #Site info stored in mem
  #Master sites list
  sites = Hash.new
  #Nested hash to hold raw data
  data = Hash.new
  vcfs.each do |vcf_file|

    name = File.basename(vcf_file, '.vcf.gz')
    data[name] = Hash.new
       
    vcfObj = Vcf.new(vcf_file)
    vcfObj.each do |key, var|

      data[name][key] = "#{var.genotypeInfo}"

      next if var.altBase=='.'
      if onlyPass
        next if var.filter != 'PASS'
        if !sites[var.uniqPos].nil?
          sites[var.uniqPos].updateQual(var.qual)
          sites[var.uniqPos].filter = var.filter
        else
          sites[var.uniqPos] = var
        end
      else
        if !sites[var.uniqPos].nil?
          sites[var.uniqPos].updateQual(var.qual)
          sites[var.uniqPos].filter = var.filter
        else
          sites[var.uniqPos] = var
        end
      end
    end
  end

  STDOUT.puts "Number of samples: #{data.length}"
  STDOUT.puts "Number of sites: #{sites.length}"

  #write header info to output

  File.open(out_file,"w") do |writer|
    #Printing VCF header information
    writer.puts "##fileformat=VCFv4.0"
    time = Time.new
    writer.puts "##fileDate=#{time.ctime}"
    writer.puts "##source=VcfPrinter v#{VERSION} r#{REVISION}"
    writer.puts "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    writer.puts "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major variant Read Depth\">"
    writer.puts "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">"
    writer.puts "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Read Depth\">"
    writer.puts "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
    writer.puts '##FORMAT=<ID=FT,Number=.,Type=String,Description="Sample Genotype Filter">'
    writer.puts '##FORMAT=<ID=AA,Number=1,Type=String,Description="Alt allele associated with variant depth">'
    writer.puts "##FILTER=<ID=low_snpqual,Description=\"SNP posterior probability lower than 0.95 (default)\">"
    writer.puts "##FILTER=<ID=low_VariantReads,Description=\"Number of variant reads is less than 3\">"
    writer.puts "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than 0.1\">"
    writer.puts "##FILTER=<ID=single_strand,Description=\"All variant reads are in a single strand direction\">"
    writer.puts "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than 6 (default)\">"
    writer.puts "##FILTER=<ID=No_data,Description=\"No valid reads on this site\">"
    writer.puts "##FILTER=<ID=No_var,Description=\"No valid variant read on this site\">"
    writer.puts "##FILTER=<ID=ReqIncl,Description=\"Site does not pass filtering requirements, but is in the list of sites required to be included\">"
    writer.puts "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t#{data.keys.join("\t")}"
         
    #Printing column names
    sites.each do |var, varObj|
      writer.puts "#{varObj.chr}\t#{varObj.coor}\t#{varObj.id}\t#{varObj.refBase}\t#{varObj.altBase}\t#{varObj.qual}\t#{varObj.filter}\t.\t#{varObj.format}"
    end
    writer.close
  end

  #Sort master sites list
  Functions.commandLineSort(out_file)


  #Loop through each sample and add columns
  data.each do |name, info|
    STDERR.puts "Processing #{name}"
    tmpFile = Tempfile.new('quickvcf')
    File.open(out_file,'r').each_line do |line|
      if line =~ /^#/
        tmpFile.puts(line)
      else
        split_line = line.split("\t")
        #vcf_line = VcfLine.new(line)
        if info.has_key?("#{split_line[0]}_#{split_line[1]}_#{split_line[3]}")
          tmpFile.puts("#{line.strip}\t#{info["#{split_line[0]}_#{split_line[1]}_#{split_line[3]}"]}")
        else
          tmpFile.puts("#{line.strip}\t.")
        end
      end
    end
    
    tmpFile.close
    FileUtils.mv(tmpFile.path, out_file)
  end
  STDERR.puts "Done!"
end


def usage
  STDERR.puts "USAGE\n\truby RunPrinter.rb \
\n\t -i \"*.vcf\" REQUIRED \n\t -o outputfile REQUIRED \n\t -l \"*.pileup\" OPTIONAL\
\n\t -p PASS (Only include PASS sites) OPTIONAL \n\t --fast OPTIONAL (Mem intensive)\
\n\t --indel OPTIONAL If input contains indels; please use this flag\
\n\t --cluster OPTIONAL \n\t -c chromosome (Option only available on cluster version) OPTIONAL\
\n\t -n number of jobs [1] (Option specific to cluster version) OPTIONAL\n"
  STDERR.puts "NOTE\n * File naming convention:\n\n\t [sample_name].vcf <.vcf file> \
\n\t [sample_name].pileup <.pileup file>"
  STDERR.puts "* Copy the settings.rb.default to settings.rb and fill in the details"
  STDERR.puts "\nCLUSTER VERSION SPECIFIC INSTRUCTIONS"
  STDERR.puts "\n\t* Fill in the cluster specific settings section in settings.rb"
  STDERR.puts "\n\t* Cluster version only accepts bgzip compressed and tabix indexed\n\t\
vcf and pileup files as input."
  STDERR.puts " \n\t* Program should be able to invoke tabix by running \n\t\
the command 'tabix' from the command line. tabix \n\t\
can be downloaded from here \n\t\
https://sourceforge.net/projects/samtools/files/tabix/"
  STDERR.puts " \n\t* How to bgzip\n\tbgzip [sample_name].vcf\n\tThis will create [sample_name].vcf.gz"
  STDERR.puts " \n\t* How to index\n\ttabix -p vcf [sample_name].vcf.gz\n\tThis will create [sample_name].vcf.gz.tbi"
  STDERR.puts " \n\t* If RunPrinter.rb was able to successfully submit jobs to cluster but you don't see \
\n\t output file, please look at the .err and .out file created by cluster command for cause\
\n\t of the error."
end

def processArgs()
  cline = "#{RUBY_EXE}  #{INSTALL_DIR}/vcfPrinter.rb "
  if ARGV.length < 4
    usage()
    Process.exit(0)
  else

    opts = GetoptLong.new(
      ["--vcfFiles","-i",GetoptLong::REQUIRED_ARGUMENT],
      ["--out","-o",GetoptLong::REQUIRED_ARGUMENT],
      ["--pass","-p",GetoptLong::NO_ARGUMENT],
      ["--indel",GetoptLong::NO_ARGUMENT],
      ["--pileups","-l",GetoptLong::OPTIONAL_ARGUMENT],
      ["--jobs","-n",GetoptLong::REQUIRED_ARGUMENT],
      ["--fast",GetoptLong::OPTIONAL_ARGUMENT],
      ["--cluster",GetoptLong::OPTIONAL_ARGUMENT],
      ["--chrom","-c",GetoptLong::REQUIRED_ARGUMENT],
      ["--usage","-h",GetoptLong::NO_ARGUMENT]
    )
    cl_args = {}
    #defaults
    cl_args['pass'] = false
    cl_args['hasPileup'] = false
    cl_args['cluster'] = false
    cl_args['fast'] = false
    cl_args['chr'] = false
    cl_args['indel'] = false
    cl_args['jobs'] = 1

    #process incoming arguments
    begin
      opts.each do |opt, arg|
        case opt
          when "--vcfFiles"
            cl_args['vcf_files'] = Dir.glob(arg)
	    cline="#{cline} -i \"#{arg}\" "
          when "--out"
            cl_args['out_file'] = arg
          when "--pass"
            cl_args['pass'] = true
            cline="#{cline} -p "
          when "--jobs"
            cl_args['jobs'] = arg.to_i
          when "--indel"
            cl_args['indel'] = true
            cline="#{cline} --indel "
          when "--chrom"
            cl_args['chr'] = arg
          when "--fast"
            cl_args['fast'] = true
          when "--cluster"
            cl_args['cluster'] = true
          when "--pileups"
	    cl_args['pileups'] = Dir.glob(arg)
            cl_args['hasPileup'] = true
	    cline="#{cline} -l \"#{arg}\" "
          when "--usage"
            usage()
            Process.exit(0)
        end
      end
      return cline, cl_args
      rescue
        usage()
        Process.exit(0)
      end
    end
  end

if __FILE__==$0
	cline, arguments = processArgs()
        arg_str = arguments.map {|k,v| "#{k}=#{v}"}.join(' , ')
        Log.info 'Start'
        Log.info "NEW RUN #{arg_str}"
        if arguments['cluster']
          cluster = VcfOnCluster.new(cline, arguments)
    	  chrjobid,chrfilename=cluster.submitClusterJobs()
	  cluster.processJobs(chrjobid, chrfilename)
          #below hack is needed bcoz chrs hash 
          #does not output chr in numerical order
          #cluster.sortVCF(arguments['out_file'])
          FileUtils.chmod(0644, arguments['out_file'])
        elsif arguments['fast']
          quickPrinter(arguments['vcf_files'], arguments['out_file'], arguments['pass'])
          if arguments['indel']
            Functions.indelProcessing(arguments['out_file'])
          else
            Functions.collapseVariants(arguments['out_file'])
          end
        else
          singleCorePrinterJob(cline, arguments['out_file'], arguments['pass'])
        end
        Log.info 'Done!'
end

