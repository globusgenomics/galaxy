#!/usr/bin/env ruby
require 'brl/util/textFileUtil'


def usage()
 $stderr.puts "extractReadNamesFromXM.rb <xm file> <out file>"
 exit(1)
end

if (ARGV.size==0 || ARGV[0]=='-h' || ARGV[1]=='--help') then
   usage()
end

inFile = ARGV[0]
outFile = ARGV[1]

r = BRL::Util::TextReader.new(inFile)
w = BRL::Util::TextWriter.new(outFile)
l = nil
r.each {|l|
   if (l=~/chr/) or (l=~/Saureus_USA300_genome/) or (l=~/gi\|87159884\|ref\|NC\_007793\.1\|/)then
   #if l=~ /^\s*(\d+)\s+(\d+\.\d+)\s+([0-9\.]+)\s+([0-9\.]+)\s+(\S+)\s+(\d+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\**?)/
      #$stderr.puts "#######"
     
     f = l.strip.split(/\s+/)
     readId = f[4]
     $stderr.puts f.size
     $stderr.puts readId
     w.puts readId
     
   end
}
r.close()
w.close()


