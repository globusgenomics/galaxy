#!/usr/bin/env ruby

require 'ccUtils.rb'
require 'brl/util/textFileUtil.rb'


module BRL; module Pash

class FastaFilter
  DEBUG =  false
  def initialize (optsHash)
		@optsHash = optsHash
		setParameters
	end
	
	def setParameters()
		@fastaFiles = @optsHash['--fastaFiles']
		@controlList = @optsHash['--readListFile']
		@outFileBaseName = @optsHash['--outFile']
		if (@optsHash['--chunk'] != nil) then
      @chunkSize = @optsHash['--chunk'].to_i
		else
      @chunkSize = 10000000
		end
		
		$stderr.puts "filtering files #{@fastaFiles} using " +
		  "the control list #{@controlList}; output base name "+
		  "is #{@outFileBaseName}, maximum chunk size is "+
		  "#{@chunkSize} sequences"
	end
		
	public
  def filterFastas()
    filterFastaFilesByReadList(@fastaFiles, @controlList, @outFileBaseName, @chunkSize)
  end
		
		
	private
	def filterFastaFilesByReadList(fastaFiles, controlList, outFileBaseName, chunkSize)
    # load the read names
    readsNamesHash = {}
    a = BRL::Util::TextReader.new(controlList)
    l = nil
    a.each {|l|
      readsNamesHash[l.strip]=1
    }
    a.close()
    fastaFileList = BRL::Pash::FileUtils.expandPatternList(fastaFiles)
    outFileIndex = 0
    outFileName = "#{outFileBaseName}.#{outFileIndex}.gz"
    outFileNextPos = 0
    outFile = BRL::Util::TextWriter.new(outFileName, "w", "gzip")
    fastaFileList.each { |fastaFileName|
      puts "processing FASTA file is #{fastaFileName}"
      fastaFile = BRL::Util::TextReader.new(fastaFileName)
      line=nil
      currentDefline = nil
      currentSequenceId = nil
      keepSequence = false
      fastaFile.each { |line|
        if (line.strip =~ /^(\s|\t)*>(\s|\t)*([^\s\t]+)/) then
          currentDefline = $~[3]
          $stderr.puts "found defline #{currentDefline}" if (DEBUG)
          if (readsNamesHash.member?(currentDefline)) then
            $stderr.puts "pos #{outFileNextPos} chunk size #{chunkSize}" if (DEBUG)
            if (outFileNextPos == chunkSize) then
              outFile.close
              outFileIndex += 1
              outFileName = "#{outFileBaseName}.#{outFileIndex}.gz"
              outFileNextPos = 0
              outFile = BRL::Util::TextWriter.new(outFileName, "w", "gzip")
            end
            outFileNextPos += 1
           #$stderr.puts "#{currentDefline} is in the control database"
            keepSequence = true
            outFile.puts "#{line}"
          else
            keepSequence = false
          end
        else # we have sequence
          if (keepSequence) then
            outFile.puts "#{line}"
          end
        end
      }
      fastaFile.close()
    }
    outFile.close()
 	end
	
	def FastaFilter.processArguments()
	# We want to add all the prop_keys as potential command line options
		optsArray =[['--fastaFiles', '-f', GetoptLong::REQUIRED_ARGUMENT],
		            ['--readListFile', '-l', GetoptLong::REQUIRED_ARGUMENT],
			          ['--chunk'     , '-c', GetoptLong::OPTIONAL_ARGUMENT],
			          ['--outFile'   , '-o', GetoptLong::REQUIRED_ARGUMENT],
			          ['--help'    , '-h', GetoptLong::NO_ARGUMENT]
		           ]
		#PROP_KEYS.each { |propName|
		#	argPropName = "--#{propName}"
		#	optsArray << [argPropName, GetoptLong::OPTIONAL_ARGUMENT]
		#}
		progOpts = GetoptLong.new(*optsArray)
		optsHash = progOpts.to_hash
		FastaFilter.usage() if(optsHash.key?('--help'));
		unless(progOpts.getMissingOptions().empty?)
			FastaFilter.usage("USAGE ERROR: some required arguments are missing") 
		end
		FastaFilter.usage() if(optsHash.empty?);
		return optsHash
	end
    
    
  def FastaFilter.usage(msg='')
			unless(msg.empty?)
				puts "\n#{msg}\n"
			end
			puts "
PROGRAM DESCRIPTION:
  Filters NCBI Fasta files by selecting only the sequences
whose trace id occur in the trace ids database.
  
COMMAND LINE ARGUMENTS:
  --fastaFiles   |  -f   => NCBI Fasta files
  --readListFile |  -l   => file containing the list of reads to be selected
  --outFile      |  -o   => output file base name
  --chunk        |  -c   => fasta file chunk size, in number of sequences
  --help         |  -h   => [optional flag] Output this usage info and exit

USAGE:
  filterNCBIFasta.rb  -f \"fasta.*\" -d NCBI.WGA.Fasta -o fasta.WGA.out -chunk 50000
";
			exit(2);
	end
end

    
end; end





#############
# Process command line options
optsHash = BRL::Pash::FastaFilter.processArguments()
# Instantiate analyzer using the program arguments
fastaFilter = BRL::Pash::FastaFilter.new(optsHash)
# Analyze this !
fastaFilter.filterFastas()
exit(0);
