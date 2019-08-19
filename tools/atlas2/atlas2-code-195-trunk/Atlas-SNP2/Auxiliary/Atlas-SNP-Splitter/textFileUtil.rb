#!/usr/bin/env ruby
# Turn on extra warnings and such

# ##############################################################################
# REQUIRED LIBRARIES
# ##############################################################################
$VERBOSE = nil
require 'zlib'
HAVE_BZ2 =  begin
              require 'bz2'
              true
            rescue LoadError => lerr
              false
            end
require 'delegate'
require 'util'
# ##############################################################################

module BRL ; module Util

	class TextReader < SimpleDelegator
		include BRL::Util
		def initialize(fileStr)
			# First, make sure we have a string obj and that the file indicated
			# exists and is readable
			unless(fileStr.kind_of?(String))
				raise(TypeError, "The file to read from must be provided as a String.")
			end
			unless(FileTest.exists?(fileStr))
				raise(IOError, "The file '#{fileStr}' doesn't exist.")
			end
			unless(FileTest.readable?(fileStr))
				raise(IOError, "The file '#{fileStr}' isn't readable.")
			end
			# Ok, figure out if we've got a gzip file or plain text file (default)
			# Once we know, use an appropriate IO delegate.
			if(Gzip.isGzippedFile?(fileStr))
				@ioObj = Zlib::GzipReader.open(fileStr)
			elsif(HAVE_BZ2 and Bzip2.isBzippedFile?(fileStr))
				@ioObj = BZ2::Reader.open(fileStr)
			else
				@ioObj = File.open(fileStr, "r")
			end
			super(@ioObj)
		end

	end # class TextReader

	class TextWriter < SimpleDelegator
		include BRL::Util
		
		GZIP_OUT, BZIP_OUT, TEXT_OUT = 0,1,2
		
		def initialize(fileStr, modeStr="w+", outputType=false)
			# First, make sure we have a string obj and that the file indicated
			# is writable
			unless(fileStr.kind_of?(String))
				raise(TypeError, "The file to write to must be provided as a String.")
			end
			if(FileTest.exists?(fileStr) and !FileTest.writable?(fileStr))
				raise(IOError, "The file '#{fileStr}' exists but isn't writable.")
			end
			# Ok, figure out if we've got a gzip file or plain text file (default)
			# Once we know, use an appropriate IO delegate.
			file = File.open(fileStr, modeStr)
			# Figure out how to zip output if asked for
			if((outputType == true) or (outputType == GZIP_OUT) or (outputType =~ /gzip/i))	# Back-compatible: true [*specifically*, not just if(outputType)] means do GZIP
				@ioObj = Zlib::GzipWriter.new(file)
			elsif(HAVE_BZ2 and ((outputType == BZIP_OUT) or (outputType =~ /bzip/i)))
				@ioObj = BZ2::Writer.new(file)
			else
				@ioObj = file
			end
			super(@ioObj)
		end
	end # class TextWriter
end ; end # module BRL ; module Util

	# ##############################################################################
	# TEST DRIVER (run this file on its own)
	# ##############################################################################
if(__FILE__ == $0)
	module TestTextFileUtil
		require 'getoptlong'

		def TestTextFileUtil.processArguments
			progOpts =
				GetoptLong.new(
					['--fileToRead', '-i', GetoptLong::REQUIRED_ARGUMENT],
					['--fileToWrite', '-o', GetoptLong::REQUIRED_ARGUMENT],
					['--writeZip', '-z', GetoptLong::NO_ARGUMENT],
					['--help', '-h', GetoptLong::NO_ARGUMENT]
				)

			optsHash = progOpts.to_hash
			return optsHash
		end

		def TestTextFileUtil.usage(msg='')
			unless(msg.empty?)
				puts "\n#{msg}\n"
			end
			puts "

	PROGRAM DESCRIPTION:

	COMMAND LINE ARGUMENTS:
	  -i		=> Location of the input file (plain or gzipped text)
	  -o    => Location of output file
	  -z    => [optional flag; default is no] Should output file be gzipped text?

	USAGE:
		textFileUtil.rb -i ./myTestFile.txt.maybeGZ.maybeNot -o ./myTestOutput.gz

	";
			exit(2);
		end # TestTextFileUtil.usage(msg='')
	end # module TestTextFileUtil

	optsHash = TestTextFileUtil.processArguments()
	if(optsHash.key?('--help') or optsHash.empty?())
		TestTextFileUtil.usage()
	end

	reader = BRL::Util::TextReader.new(optsHash['--fileToRead'])
	writer = BRL::Util::TextWriter.new(optsHash['--fileToWrite'], optsHash.key?('--writeZip'))

	writer.write("OUTPUT OF TextReader & TextWriter TEST DIRVER in textFileUtil.rb\n")

	reader.each {
		|line|
		writer.write(line)
	}

	reader.close()
	writer.close()
end # if(__FILE__ == $0)
