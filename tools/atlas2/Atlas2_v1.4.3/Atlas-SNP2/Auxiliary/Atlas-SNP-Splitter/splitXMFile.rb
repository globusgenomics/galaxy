#!/usr/bin/env ruby
require 'textFileUtil'
require 'util'

	
class SplitXMOutput
  DEBUG = false
  DEBUG_HASH = false
	DEBUG_GAP = false
	DEBUG_SPLIT = false
	DEBUG_BS  = false
  def initialize(optsHash)
		@optsHash = optsHash
		setParameters()
	end  
  
  def setParameters()
		@xmFile = @optsHash['--xmFile']
		if (!@optsHash.key?('--splitByChrom') && !@optsHash.key?('--numberOfParts')) then
			$stderr.puts "You need to specify a splitting method.\n"
			usage()
			exit(2)
		end
		if (@optsHash.key?('--splitByChrom') && @optsHash.key?('--numberOfParts')) then
			$stderr.puts "You need to specify exactly one splitting method.\n"
			usage()
			exit(2)
		end
		@splitByChromFlag = true
		@numberOfParts = 0
		if (@optsHash.key?('--numberOfParts')) then
			@splitByChromFlag = false
			@numberOfParts = @optsHash['--numberOfParts'].to_i
			if (@numberOfParts<=1) then
				$stderr.puts "The number of parts should be an integer greater than one."
				usage()
				exit(2)
			end
		end
		@outputFileRoot = @optsHash['--outputFileRoot']
  end
  
  def work()
		if (@splitByChromFlag) then
			splitByChrom()
		else
			splitInEqualSizedNumberOfParts()
		end
  end
  
  
  def splitByChrom()
		#initialize file handles hash
		outFileHandlesHash = {}
  
		l = nil
		@currentFileHandle = nil
		xmReader = File.open(@xmFile, "r")
		xmReader.each {|l|
			f = l.strip.split(/\s+/)
			$stderr.puts f.join(";\t") if (DEBUG)
			
			if (f.size>6) then
				if (f[f.size-1]=="*") then
					chromosomeName = f[f.size-5]
				else
					chromosomeName = f[f.size-4]
				end
				
				if (chromosomeName =~ /(.*)_(\d+)_(\d+)/) then
          chromosomeName = $1
				end
			      
				#if (chromosomeName !~ /chr/) then
				#	puts "problem line: #{l}\n#{f.join(";")}"
				#	exit(2)
				#end
				$stderr.puts "def line on chrom #{chromosomeName}" if (DEBUG)
				if (outFileHandlesHash.key?(chromosomeName)) then
					@currentFileHandle = outFileHandlesHash[chromosomeName]
				else
					@currentFileHandle = File.open("#{@outputFileRoot}.#{chromosomeName}", "w")
					outFileHandlesHash[chromosomeName] = @currentFileHandle
				end
				@currentFileHandle.print l
				# summary line
			else
				# point changes line
				@currentFileHandle.print l
			end
		}
		xmReader.close()
  
		k = nil
		outFileHandlesHash.keys.each {|k|
			outFileHandlesHash[k].close()
		}
  end
  
  
  
  def checkSameChromosome(xmFile)
    $stderr.puts "check same chromo"
		currentChrom = nil
		chromosomeName = nil
		l = nil
		@currentFileHandle = nil
		xmReader = File.open(xmFile, "r")
		xmReader.each {|l|
                  #$stderr.puts "###"
			f = l.strip.split(/\s+/)
			$stderr.puts f.join(";\t") if (DEBUG)
			
			if (f.size>6) then
				if (f[f.size-1]=="*") then
					chromosomeName = f[f.size-5]
				else
					chromosomeName = f[f.size-4]
				end
				
				if (chromosomeName =~ /(.*)_(\d+)_(\d+)/) then
          chromosomeName = $1
				end
				
				$stderr.puts "def line on chrom #{chromosomeName}" if (DEBUG)
				if (chromosomeName != currentChrom) then
					if (currentChrom==nil) then
						currentChrom = chromosomeName
						$stderr.puts currentChrom
					else
						$stderr.puts "#{xmFile} contains mappings onto chromomosomes #{currentChrom} and #{chromosomeName}"
						exit(2)
					end
				end
			end
		}
		xmReader.close()	
		return currentChrom
	end
  
  
  def buildReadMappingHash(xmFile)
		currentChrom = nil
		chromosomeName = nil
		l = nil
	
		xmReader = File.open(xmFile, "r")
		startPos = nil
		stopPos = nil
		readId = nil
		currentReadInfo = nil
		@readMappingHash = {}
		@readMappingArray = []
		xmReader.each {|l|
			f = l.strip.split(/\s+/)
			if (f.size>6) then
				$stderr.puts f.join(";\t") if (DEBUG_HASH)
				readId = f[4]
				
				
				if (f[f.size-1]=="*") then
					chromosomeName = f[f.size-5]
          if (f[f.size-6]=="C") then
            startPos = f[f.size-2].to_i
            stopPos = f[f.size-3].to_i
          else
            startPos = f[f.size-4].to_i
            stopPos = f[f.size-3].to_i
          end
				else
          chromosomeName = f[f.size-4]
          if (f[f.size-5]=="C") then
            startPos = f[f.size-1].to_i
            stopPos = f[f.size-2].to_i
          else
            startPos = f[f.size-3].to_i
            stopPos = f[f.size-2].to_i
          end
				end
				
				if (chromosomeName =~ /(.*)_(\d+)_(\d+)/) then
          chromosomeName = $1
          startPos += $2.to_i
          stopPos += $2.to_i
				end
				
				$stderr.puts "def line for #{readId} on chrom #{chromosomeName}: #{startPos} #{stopPos}" if (DEBUG_HASH)
				currentReadInfo = @readInfoStruct.new(readId, startPos, stopPos)
				@readMappingHash[readId] = currentReadInfo
				@readMappingArray.push(currentReadInfo)
			end
		}
		xmReader.close() 
  end
  
  def splitInEqualSizedNumberOfParts()
		# verify all reads are on the same chromosome
		currentChrom = checkSameChromosome(@xmFile)
		# extract read info
		@readInfoStruct = Struct.new("ReadInfo", :readId, :startPos, :stopPos)
		$stderr.puts "building read mappings array #{Time.now()}"
		buildReadMappingHash (@xmFile)
		# sort it by chrom start, stop
		$stderr.puts "Sorting the read array \n#{Time.now()}"
		@readMappingArray.sort! { |a,b| a.startPos <=> b.startPos}
		$stderr.puts "Finished sorting the read array \n#{Time.now()}"
		# use min pos, max pos, to guide split in multiple parts
		quasiEqualSplittingByNumberOfReads()
		# break it there, start populating new list
		splitByPartNumber()
	end
	
	
	
	def splitByPartNumber
		outFileHandlesHash = {}
		xmReader = File.open(@xmFile, "r")
		startPos = nil
		stopPos = nil
		readId = nil
		currentReadInfo = nil
		@currentFileHandle = nil
		outFileIndex = 0
		
		
		minPart = nil
		maxPart = nil
		found = nil
		mid = nil
		
		
		xmReader.each {|l|
			f = l.strip.split(/\s+/)
			if (f.size>6) then
				$stderr.puts f.join(";\t") if (DEBUG_SPLIT)
				chromosomeName = f[f.size-4]
				readId = f[4]
				if (!@readMappingHash.key?(readId)) then
					$stderr.puts "read #{readId} not found in hash"
					exit(2)
				end
				readInfo = @readMappingHash[readId]
				
				
				# binary search on parts
				minPart = 0
				maxPart = @stopPart.size-1
				
				found = false
				while ( !found && (minPart <= maxPart)) do
					mid = (minPart+maxPart)/2
					$stderr.puts "m=#{minPart} M=#{maxPart} mid=#{mid} compare #{readInfo.startPos} with #{@startPart[mid]} " if (DEBUG_BS)
					if (readInfo.startPos <@startPart[mid]) then
						maxPart = mid-1
					elsif (readInfo.startPos > @stopPart[mid]) then
						minPart = mid+1
					else
						found = true
						outFileIndex = mid
					end
				end
				
				if (!found) then
					$stderr.puts "Part not found for read #{readInfo}, setting it to part 0 by default "
					outFileIndex = 0
				end
				
				$stderr.puts "def line for #{readId} on chrom #{chromosomeName}: #{readInfo.startPos} #{readInfo.stopPos} ends up on part #{outFileIndex}"  if (DEBUG_HASH)
				
				if (!outFileHandlesHash.key?(outFileIndex)) then
					outFileHandlesHash[outFileIndex] = File.open("#{@outputFileRoot}.part.#{outFileIndex}", "w")
				end
				@currentFileHandle = outFileHandlesHash[outFileIndex]
				@currentFileHandle.print l
			else
				@currentFileHandle.print l
			end
		}
		xmReader.close() 	
		
		k = nil
		outFileHandlesHash.keys.each {|k|
			outFileHandlesHash[k].close()
		}	
	end
	
	
	def quasiEqualSplitting ()
	  # traverse sort array until reaching split point - 10%
		# traverse up until finding a gap !!!
		minPoint = @readMappingArray[0].startPos
		maxPoint = @readMappingArray[@readMappingArray.size-1].stopPos+1000
		# setup gap targets
		$stderr.puts "splitting mappings between #{minPoint} and #{maxPoint} in #{@numberOfParts} parts"
		
		i = nil
		startGap = []
		stopGap = []
		targetGap = []
		bestGapPos = nil
		bestGapDistance = nil
		partSize = (maxPoint - minPoint) / @numberOfParts
		actualGaps = []
		1.upto(@numberOfParts-1) {|i|
			startGap.push(minPoint + partSize*i - partSize/3)
			stopGap.push (minPoint + partSize*i + partSize/3)
			targetGap.push(minPoint + partSize*i)
		}
		
		currentMaxStopPos = minPoint
		readStruct = nil
		
		currentTargetGap = targetGap[0]
		currentGapStart = startGap[0]
		currentGapStop = stopGap[0]
		currentGapIndex = 0
		bestGapPos = nil
		bestGapDistance = 1000000000
		gapDistance = nil
		gapPos = nil
		
		@readMappingArray.each {|readStruct|
			if (readStruct.startPos > currentGapStop) then
				actualGaps.push(bestGapPos)
				currentGapIndex += 1
				if (currentGapIndex >= @numberOfParts-1) then
					break
				else
					currentTargetGap = targetGap[currentGapIndex]
					currentGapStart = startGap[currentGapIndex]
					currentGapStop = stopGap[currentGapIndex]
					bestGapPos = nil
					bestGapDistance = 1000000000			
				end
			end
			if (readStruct.startPos > currentMaxStopPos) then
				# gap found
				if (readStruct.startPos >= currentGapStart && readStruct.startPos <= currentGapStop) then
					gapPos = (readStruct.startPos+currentMaxStopPos)/2
					gapDistance = currentTargetGap - gapPos
					if (gapDistance < 0) then
						gapDistance = - gapDistance
					end
					$stderr.puts "[#{currentGapStart}\t#{currentTargetGap}\t#{currentGapStop}] gap pos #{gapPos} dist #{gapDistance}" if (DEBUG_GAP)
					if (gapDistance < bestGapDistance) then
						# replace old best gap
						bestGapDistance = gapDistance
						bestGapPos = gapPos
						$stderr.puts "updated best gap to #{bestGapPos}\t#{bestGapDistance}" if (DEBUG_GAP)
					end
				end
			end
			if (readStruct.stopPos > currentMaxStopPos) then
				currentMaxStopPos = readStruct.stopPos
			end
		}
		
		# did we process the last gap ?
		if (currentGapIndex<@numberOfParts) then
			actualGaps.push(bestGapPos)
		end
		
		# finalize analysis; if a gap not found, reduce number of parts
		@startPart = [minPoint]
		@stopPart = []
		
		gapPos = nil
		actualGaps.each {|gapPos|
			if (gapPos !=nil) then
				@stopPart.push (gapPos)
				@startPart.push(gapPos)
			end
		}
		@stopPart.push(maxPoint)
		
		0.upto(@startPart.size-1) {|i|
			$stderr.puts "interval #{i}: [#{@startPart[i]}, #{@stopPart[i]}]" 	
		}
	end
	

	def quasiEqualSplittingByNumberOfReads ()
	  
	  numberOfReadsInCurrentBin = 0
	  targetNumberOfReadsInBin = @readMappingArray.size / @numberOfParts+1
	  $stderr.puts "Splitting input in #{@numberOfParts} parts; the target number of reads per split output part #{targetNumberOfReadsInBin}"
	  currentMaxStopPos = @readMappingArray[0].stopPos
		actualGaps = []
		
		@readMappingArray.each {|readStruct|
			numberOfReadsInCurrentBin +=1	
			if (readStruct.startPos > currentMaxStopPos) then
				# gap found
				if (numberOfReadsInCurrentBin>= targetNumberOfReadsInBin) then
					$stderr.puts"Current bin has #{numberOfReadsInCurrentBin}"
					actualGaps.push( (readStruct.startPos + currentMaxStopPos)/2)
					numberOfReadsInCurrentBin = 0
				end
			end
			if (readStruct.stopPos > currentMaxStopPos) then
				currentMaxStopPos = readStruct.stopPos
			end
		}
				
		# finalize analysis; if a gap not found, reduce number of parts
		@startPart = [ @readMappingArray[0].startPos]
		@stopPart = []
		
		gapPos = nil
		actualGaps.each {|gapPos|
			@stopPart.push (gapPos)
			@startPart.push(gapPos)
		}
		@stopPart.push(@readMappingArray[@readMappingArray.size-1].stopPos+3000)
		
		0.upto(@startPart.size-1) {|i|
			puts "interval #{i}: [#{@startPart[i]}, #{@stopPart[i]}]" 	
		}
	end

	
  def SplitXMOutput.processArguments()
		# We want to add all the prop_keys as potential command line options
		optsArray =	[ ['--xmFile',     							'-x', GetoptLong::REQUIRED_ARGUMENT],
									['--outputFileRoot',					'-o', GetoptLong::OPTIONAL_ARGUMENT],
									['--splitByChrom',					  '-C', GetoptLong::OPTIONAL_ARGUMENT],
									['--numberOfParts',						'-n', GetoptLong::OPTIONAL_ARGUMENT],
									['--help',           					'-h', GetoptLong::NO_ARGUMENT]
								]
		
		progOpts = GetoptLong.new(*optsArray)
		optsHash = progOpts.to_hash
		SplitXMOutput.usage() if(optsHash.key?('--help'));
		
		unless(progOpts.getMissingOptions().empty?)
			SplitXMOutput.usage("USAGE ERROR: some required arguments are missing") 
		end
	
		SplitXMOutput.usage() if(optsHash.empty?);
		return optsHash
	end
	
	def SplitXMOutput.usage(msg='')
		unless(msg.empty?)
			puts "\n#{msg}\n"
		end
		puts "
PROGRAM DESCRIPTION:
    Splits a cross_match output file, either by chromosomes,
  or in read sets mapping on equal size chromosome regions.
  The recommended usage is to first perform the split by chromosome,
  then split the output for each chromosome in equal size parts.
  NOTE: to perform the second type of split, the cross_match file should
  contain mappings to one chromosome only.
  
COMMAND LINE ARGUMENTS:
  --xmFile                 | -x   => cross_match file
  --splitByChrom           | -C   => split input by chromosome
  --numberOfParts          | -n   => split input in a number of parts with
                                     mappings on genomic regions of the same size
  --outputFileRoot         | -o   => output file root                                    
  --help           |-h   => [optional flag] Output this usage info and exit

EXAMPLES:
  splitXMFile.rb  -x xmOutput -o xmByChrom -C
  splitXMFile.rb -x xmByChrom.chr1 -n 10 -o xmByChrom.chr1.10parts 
";
			exit(2);
	end
end


########################################################################################
# MAIN
########################################################################################

# Process command line options
optsHash = SplitXMOutput.processArguments()
# Instantiate analyzer using the program arguments
xmSplitter = SplitXMOutput.new(optsHash)
# Analyze this !
xmSplitter.work()
exit(0);
