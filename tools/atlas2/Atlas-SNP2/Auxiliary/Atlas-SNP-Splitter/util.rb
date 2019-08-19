#!/usr/bin/env ruby

# Turn off extra warnings and such
$VERBOSE = nil

require 'getoptlong'
require 'fcntl'
require 'zlib'
require 'cgi'

module BRL ; module Util
  FATAL, OK, OK_WITH_ERRORS, FAILED, USAGE_ERR = 1,0,10,20,16
  BLANK_RE = /^\s*$/
  COMMENT_RE = /^\s*#/

  def setLowPriority()
    begin
      Process.setpriority(Process::PRIO_USER, 0, 19)
    rescue
    end
    return
  end

  class Job
    def Job::makeJobTicket()
      pid = Process.pid
      time = Time::now.to_i
      return "#{time}-#{pid}"
    end
  end

  class MemoryInfo
    @@procStatusFH = nil

    def MemoryInfo.getMemUsageStr()
      memStr = ''
      if(@@procStatusFH.nil?)
        begin
          @@procStatusFH = File.open("/proc/#{$$}/status")
        rescue
          return nil
        end
      end
      @@procStatusFH.rewind
      @@procStatusFH.readlines.each { |line|
        line.strip!
        next unless(line =~ /^VmSize:\s+(.*)$/)
        memStr = $1
        break
      }
      return memStr.strip
    end

    def MemoryInfo.getMemUsagekB()
      memStr = MemoryInfo.getMemUsageStr()
      return -1 if(memStr.nil? or memStr.empty?)
      memStr =~ /^(\d+)\s+(\S+)$/
      memUsage, units = $1.to_i, $2.strip.downcase
      if(units =~ /kb/)
        return memUsage
      elsif(units =~ /mb/)
        return memUsage * 1000
      elsif(units =~ /gb/)
        return memUsage * 1000 * 1000
      else # bytes? wtf
        return -1
      end
    end
  end

  module Bzip2
    BZIP2_EXE = 'bzip2 -t '

    def Bzip2.isBzippedFile?(fileStr)
      fullFilePath = File.expand_path(fileStr)
      return false if((!fullFilePath.kind_of?(String)) or !FileTest.exists?(fullFilePath) or !FileTest.readable?(fullFilePath))
      fullFilePath.gsub!(/\"/, '\\"')
      cmdStr1 = "#{BRL::Util::Bzip2::BZIP2_EXE} \"#{fullFilePath}\" > /dev/null 2> /dev/null"
      begin
        system(cmdStr1)
        exitCode = $?
        return (exitCode == 0) ? true : false
      rescue
        raise "\n\nERROR: can't execute sub-processes??? (BRL::Util.Util#Bzip.isBzippedFile?)\n\n"
      end
    end

  end  # module Bzip2

  module Gzip

    GZIP_EXE = 'gzip -t '
    GZIP_INFO = 'gzip -l'
    GZIP_INFO_RE = /\s*compressed\s*uncompressed\s*ratio\s*uncomp\S+\s*\r?\n\s*(\d+)\s*(\d+)\s*([0-9\.]+)\%\s*(\S+)/

    def Gzip.isGzippedFile?(fileStr)
      fullFilePath = File.expand_path(fileStr)
      return false if((!fullFilePath.kind_of?(String)) or !FileTest.exists?(fullFilePath) or !FileTest.readable?(fullFilePath))
      fullFilePath.gsub!(/\"/, '\\"')
      cmdStr1 = "#{BRL::Util::Gzip::GZIP_EXE} \"#{fullFilePath}\" > /dev/null 2> /dev/null "
      begin
        system(cmdStr1)
        exitCode = $?
        return (exitCode == 0) ? true : false
      rescue Exception => err
        raise "\n\nERROR: can't execute sub-process '#{BRL::Util::Gzip::GZIP_EXE}' (BRL::Util.Util#Gzip.isGzippedFile?)\nBACKTRACE:#{err.message}\n" + err.backtrace.join("\n") + "\n\n"
      end
    end

    def Gzip.getInfo(fileStr)
      fullFilePath = File.expand_path(fileStr)
      return nil if((!fileStr.kind_of?(String)) or !FileTest.exists?(fullFilePath) or !FileTest.readable?(fullFilePath))
      fullFilePath.gsub!(/\"/, '\\"')
      cmdStr = "#{BRL::Util::Gzip::GZIP_INFO} \"#{fullFilePath}\" "
      begin
        cmdOutput = `#{cmdStr}`
      rescue => err
        raise(SystemCallError, "\nERROR: 'gzip -l' doesn't work! do you have it installed and reachable from sh? Is your gzip up to date?\nThis command returned:]\n  #{cmdOutput}\nError details: #{err.message}\n" + err.backtrace.join("\n")) if(cmdOutput.nil? or cmdOutput.empty?)
        raise(SystemCallError, "\nERROR: '#{fileStr}' is not a valid gzipped file!") if(cmdOutput =~ /not in gzip format/)
      end
      if(cmdOutput =~ BRL::Util::Gzip::GZIP_INFO_RE)
        return [$1, $2, $3.to_f/100.0, $4]
      else
        raise(SystemCallError, "\nERROR: 'gzip -l' doesn't work! do you have it installed and reachable from sh? Is your gzip up to date?\nThis command returned:\n>>\n#{cmdOutput}\n<<")
      end
    end

    def Gzip.empty?(fileStr)
      fullFilePath = File.expand_path(fileStr)
      if((!fileStr.kind_of?(String)) or !FileTest.exists?(fullFilePath) or !FileTest.readable?(fullFilePath))
        raise(SystemCallError, "\nERROR: the file isn't the name of a gzippe file or the file doesn't exist or you don't have permission to read it!\n    Bad file: #{fileStr}\nBACKTRACE:#{err.message}\n" + err.backtrace.join("\n") + "")
      end
      testByte = nil
      begin
        # Try to read just 1 byte from gzip file
        # If fails, that's because the compressed file was empty!
        gzIo = Zlib::GzipReader.open(fullFilePath)
        testByte = gzIo.readchar
      rescue
      end
      return testByte.nil? # yes, compressed file was empty!
    end

    def Gzip.getCompressedSize(fileStr)
      return BRL::Util::Gzip::getGzipInfo(fileStr)[0]
    end

    def Gzip.getUncompressedSize(fileStr)
      return BRL::Util::Gzip::getGzipInfo(fileStr)[1]
    end

    def Gzip.getRatio(fileStr)
      return BRL::Util::Gzip::getGzipInfo(fileStr)[2]
    end

    def Gzip.getFileName(fileStr)
      return BRL::Util::Gzip::getGzipInfo(fileStr)[3]
    end
  end # module Gzip
end ; end # module BRL ; module Util

# ##############################################################################
# ADD TO EXISTING CLASSES/MODULES
# ##############################################################################

class Object
  def eputs(*args)
    eprint(*args)
    $stderr.print "\n"
    return
  end

  def eprint(*args)
    args.each { |arg| $stderr.print arg }
    return
  end

  def dputs(*args)
    if($DEBUG)
      eputs(*args)
    end
    return
  end

  def dprint(*args)
    if($DEBUG)
      eprint(*args)
    end
    return
  end

  def initAttributes(attrList)
    attrList.each { |attrSym|
      tmpSelfClass = class << self
                       self
                     end
      tmpSelfClass.module_eval { attr_accessor(attrSym) }
      self.send("#{attrSym}=", nil)
    }
    return
  end

  ########################################################################
  # * *Function*: A deep copy mechanism.
  #
  # * *Usage*   : <tt>  deep_clone( [ ["a"], ["b", "c"]] )  </tt>
  # * *Args*    :
  #   - +Object+ -> The object to clone
  # * *Returns* :
  #   - +Object+ -> A clone of the given object.
  # * *Throws* :
  #   - +none+
  ########################################################################
  def deep_clone
      Marshal::load(Marshal.dump(self))
  end
end

class Symbol
  def <=>(otherSym)
    (retVal = (self.to_s.downcase <=> otherSym.to_s.downcase))
    (retVal = (self.to_s <=> otherSym.to_s)) if(retVal == 0)
    retVal
  end
end

# We'd like the options available in some Hash object we can query at any time
class GetoptLong
  def to_hash
    return @asHash if(self.instance_variables.include?("@asHash") and !@asHash.nil?)
    saveSettings()
    retHash = {};
    self.each { |optName, optValue|
      retHash[optName] = optValue
    }
    # Do url-unescaping if it looks like the value has at least one actual %-style url encoded value.
    retHash.each_key { |kk| retHash[kk] = CGI.unescape(retHash[kk]) if(retHash[kk] =~ /%[0-9a-fA-F][0-9a-fA-F]/) ; }
    restoreSettings()
    @asHash = retHash
    return retHash
  end

  def getMissingOptions(butRequired=true)
    retVal = []
    self.to_hash()
    self.getCanonicalNames().each { |optName|
      if(butRequired)
        retVal << optName if(!@asHash.key?(optName) and @argument_flags[optName] == REQUIRED_ARGUMENT)
      else # any missing
        retVal << optName unless(@asHash.key?(optName))
      end
    }
    return retVal
  end

  def getCanonicalNames()
    retVal = {}
    @canonical_names.each_value { |optName|  retVal[optName] = nil }
    return retVal.keys
  end

  def saveSettings()
    @old_canonical_names = @canonical_names unless(@canonical_names.nil?)
    @old_argument_flags = @argument_flags unless(@argument_flags.nil?)
    @old_non_option_arguments = @non_option_arguments unless(@non_option_arguments.nil?)
    @old_rest_singles = @rest_singles unless(@rest_singles.nil?)
    return
  end

  def restoreSettings()
    @canonical_names = @old_canonical_names unless(@old_canonical_names.nil?)
    @argument_flags = @old_argument_flags unless(@old_argument_flags.nil?)
    @non_option_arguments = @old_non_option_arguments unless(@old_non_option_arguments.nil?)
    @rest_singles = @old_rest_singles unless(@old_rest_singles.nil?)
    return
  end

end # class GetoptLong

# Add overlap detection and length/size to Range object
class Range
  def size
    self.last - self.first + (self.exclude_end?() ? 0 : 1)
  end

  def exclude_max?()
    if(exclude_end?())
      if(self.first < self.last)
        return true
      end
    end
    return false
  end

  def exclude_min?()
    if(exclude_end?())
      if(self.first > self.last)
        return true
      end
    end
    return false
  end

  def <=>(aRange)
    selfMin = (self.first < self.last ? self.first : self.last)
    arMin = (aRange.first < aRange.last ? aRange.first : aRange.last)
    selfMax = (self.first > self.last ? self.first : self.last)
    arMax = (aRange.first > aRange.last ? aRange.first : aRange.last)
    selfExclMax, arExclMax = self.exclude_max?(), aRange.exclude_max?()
    selfExclMin, arExclMin = self.exclude_min?(), aRange.exclude_min?()

    # Compare min sentinals
    retVal = selfMin <=> arMin
    # Resolve ties where exclude_mins differs
    if( retVal == 0 )
      if( selfExclMin and !arExclMin )
        retVal = 1
      elsif( !selfExclMin and arExclMin )
        retVal = -1
      end
    end
    # If *still* a tie, then min sentinals are equal and exclude_mins are same
    if( retVal == 0 )
      retVal = selfMax <=> arMax
      # Resolve ties where exlcude_max differs
      if(retVal == 0 )
        if( selfExclMax and !arExclMax )
          retVal = -1
        elsif( !selfExclMax and arExclMax )
          retVal = 1
        end
      end
    end
    return retVal
  end

  def contains?(aRange)
    raise(TypeError, "\nERROR: containsRange?() requires Range-like methods to be available: first(), last(), ===, exclude_end?") unless(aRange.respond_to?("===") and aRange.respond_to?("first") and aRange.respond_to?("last") and aRange.respond_to?("exclude_end?"))
    selfMin, arMin = (self.first < self.last ? self.first : self.last), (aRange.first < aRange.last ? aRange.first : aRange.last)
    selfMax, arMax = (self.first > self.last ? self.first : self.last), (aRange.first > aRange.last ? aRange.first : aRange.last)

    # Check smaller (min) end of aRange
    if( self.exclude_min?() )
      if( aRange.exclude_min?() )
        return false unless(arMin >= selfMin)
      else
        return false unless(arMin > selfMin)
      end
    else
      return false unless(arMin >= selfMin)
    end
    # Check larger (max) end of aRange
    if( self.exclude_max?() )
      if( aRange.exclude_max?() )
        return false unless(arMax <= selfMax)
      else
        return false unless(arMax < selfMax)
      end
    else
      return false unless(arMax <= selfMax)
    end
    # Must contain aRange
    return true
  end

  def within?(aRange)
    raise(TypeError, "\nERROR: argument must respond to contains?() and other range-like methods.") unless(aRange.respond_to?("===") and aRange.respond_to?("first") and aRange.respond_to?("last") and aRange.respond_to?("exclude_end?") and aRange.respond_to?('contains?') )
    return aRange.contains?(self)
  end

  def overlaps?(aRange)
    raise(TypeError, "\nERROR:rangesOverlap?() requires Range-like methods to be available: first(), last(), ===, exclude_end?") unless(aRange.respond_to?("===") and aRange.respond_to?("first") and aRange.respond_to?("last") and aRange.respond_to?("exclude_end?"))
    selfMin, arMin = (self.first < self.last ? self.first : self.last), (aRange.first < aRange.last ? aRange.first : aRange.last)
    selfMax, arMax = (self.first > self.last ? self.first : self.last), (aRange.first > aRange.last ? aRange.first : aRange.last)

    # Check smaller (min) end of aRange
    if( self.exclude_max?() )
      if( aRange.exclude_min?() )
        return false unless(arMin <= selfMax)
      else
        return false unless(arMin < selfMax)
      end
    else
      return false unless(arMin < selfMax)
    end
    # Check smaller (min) of self
    if( aRange.exclude_max?() )
      if( self.exclude_min?() )
        return false unless(selfMin <= arMax)
      else
        return false unless(selfMin < arMax)
      end
    else
      return false unless(selfMin < arMax)
    end
    # Must overlap
    return true
  end

  # Only works for methods seen so far
  alias_method :length, :size
  alias_method :containsRange?, :contains?
  alias_method :rangesOverlap?, :overlaps?
end # class Range

class Numeric
  def commify()
    return self.to_s.commify
  end
end

# Add NaN class variable to Float
class Float
  NaN = 0.0 / 0.0
end

# Add machine-specific native Integer sizes (before Bignum gets involved)
class Fixnum
  MAX32 = 2**31 - 1
  MIN32 = -MAX32 - 1
  MAX64 = 2**63 - 1
  MIN64 = -MAX64 - 1
  N_BYTES_NATIVE =[42].pack('i').size
  N_BITS_NATIVE = N_BYTES_NATIVE * 8
  MAX_NATIVE = 2 ** (N_BITS_NATIVE - 1) - 1
  MIN_NATIVE = -MAX_NATIVE - 1
end

class Integer
  MAX32 = Fixnum::MAX32
  MIN32 = Fixnum::MIN32
  MAX64 = Fixnum::MAX64
  MIN64 = Fixnum::MIN64
  N_BYTES_NATIVE = Fixnum::N_BYTES_NATIVE
  N_BITS_NATIVE = Fixnum::N_BITS_NATIVE
  MAX_NATIVE = Fixnum::MAX_NATIVE
  MIN_NATIVE = Fixnum::MIN_NATIVE
end

# Add safeMkdir to Dir class
class Dir
  def Dir.safeMkdir(dirToMake)
    Dir.mkdir(dirToMake) unless(File.exist?(dirToMake))
  end

  def Dir.recursiveSafeMkdir(dirToMake)
    dirNodes = dirToMake.split('/')
    currDir = ''
    dirNodes.each {
      |dirNode|
      currDir << "#{dirNode}/"
      Dir.mkdir(currDir) unless(File.exist?(currDir))
    }
  end
end

class File
  # NFS-Safe flock()
  # By Matz, with small robustification by ARJ
  alias __flock flock
  def flock(cmdFlags)
    # We trying to block or what?
    icmd = ((cmdFlags & LOCK_NB) == LOCK_NB) ? Fcntl::F_SETLK : Fcntl::F_SETLKW
    type =
      case(cmdFlags & ~LOCK_NB)
        when LOCK_SH
          Fcntl::F_RDLCK
        when LOCK_EX
          Fcntl::F_WRLCK
        when LOCK_UN
          Fcntl::F_UNLCK
        else
          raise ArgumentError, cmd.to_s
      end
    flock = [type, 0, 0, 0, 0].pack("ssqqi")
    begin
      fcntl(icmd, flock)
    rescue Errno::EBADF => err
      raise ArgumentError, "this file is not open for writing. Cannot LOCK_EX files open for reading only.\n\n"
    end
  end

  def getLock()
    loop {
      begin
        self.flock(File::LOCK_EX | File::LOCK_NB)
        break
      rescue Errno::EAGAIN => err
        # This is ok, the resourse is in use
        sleep(1)
      end
    }
  end
end

class Time
  def to_rfc822()
    if(Time.instance_methods.include?('rfc822'))
      return self.rfc822()
    else
      return self.strftime("%a, %d %b %Y %H:%M:%S ") + sprintf("%+03i00", (Time.now.utc_offset / 60 / 60))
    end
  end
end

class String
  def commify()
    return self.gsub(/,/, '').gsub(/(\d)(?=\d{3}+(?:\.|$))(\d{3}\..*)?/, '\1,\2')
  end

  ########################################################################
  # * *Function*: A fuzzy matching mechanism (still a bit rough)
  #
  # * *Usage*   : <tt>  "Alexsander".fuzzy_match( "Aleksander" )  </tt>
  # * *Args*    :
  #   - +str_in+ -> The string to compare against
  # * *Returns* :
  #   - +score+ -> A score from 0-1, based on the number of shared edges (matches of a sequence of characters of length 2 or more)
  # * *Throws* :
  #   - +none+
  ########################################################################
  def fuzzy_match( str_in )
  # The way this works:
  #   Converts each string into a "graph like" object, with edges
  #       "alexsander" - > [ alexsander, alexsand, alexsan ... lexsand ... san ... an, etc ]
  #       "aleksander" - > [ aleksander, aleksand ... etc. ]
  #   Perform match, then remove any subsets from this matched set (i.e. a hit on "san" is a subset of a hit on "sander")
  #       Above example, once reduced -> [ ale, sander ]
  #   See's how many of the matches remain, and calculates a score based
  #   On how many matches, their length, and compare to the length of the larger of the two words
    return 0 if str_in == nil
    return 1 if self == str_in
    # Make a graph of each word (okay, so its not a true graph, but is similar
    graph_A = Array.new
    graph_B = Array.new

    # "graph" self
    last = self.length
    (0..last).each{ |ff|
        loc  = self.length
        break if ff == last - 1
        wordB = (1..(last-1)).to_a.reverse!
        wordB.each{ |ss|
            break if ss == ff
            graph_A.push( "#{self[ff..ss]}" )
        }
    }

    # "graph" input string
    last = str_in.length
    (0..last).each{ |ff|
        loc  = str_in.length
        break if ff == last - 1
        wordB = (1..(last-1)).to_a.reverse!
        wordB.each{ |ss|
            break if ss == ff
            graph_B.push( "#{str_in[ff..ss]}" )
        }
    }

    # count how many of these "graph edges" we have that are the same
    matches = Array.new
    graph_A.each{ |aa|
        matches.push( aa ) if( graph_B.include?( aa ) )
    }

    matches.sort!{ |x,y| x.length <=> y.length }  # For eliminating subsets, we want to start with the smallest hits

    # eliminate any subsets
    mclone = matches.dup
    mclone.each_index { |ii|
        reg = Regexp.compile( mclone[ii] )
        count = 0.0
        matches.each{ |xx|
            count += 1 if xx =~ reg
        }
        matches.delete(mclone[ii]) if count > 1
    }

    score = 0.0
    matches.each{ |mm| score += mm.length }

    self.length > str_in.length ? largest = self.length : largest = str_in.length
    return score/largest
  end
end

class Array
  def lastIndex()
    return nil if self.empty?
    return (self.size-1)
  end
end # class Array

class CGI
  class << self
    alias :escape_orig :escape
    def escape(string)
      retVal = CGI::escape_orig(string.to_s)
      retVal.gsub!(/\+/, '%20')
      return retVal
    end
  end
end
