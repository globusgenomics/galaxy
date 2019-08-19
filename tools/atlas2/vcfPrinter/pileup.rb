require 'pileup_variant.rb'


class Pileup < Hash

        attr_accessor :pileup

        def initialize(pileup_file, *args)
          if args.length != 0 
            groups = %r/(\w+):(\d+)-(\d+)/.match(args[0])
            chr = groups[1]
            start_pos = groups[2]
            end_pos = groups[3]
          end

          if File.fnmatch("*.gz", pileup_file)
            output = `tabix #{pileup_file} #{chr}:#{start_pos.to_i}-#{end_pos.to_i}`.strip.split("\n")
          else
            output = File.open(pileup_file,'r').readlines
          end

          output.each do |out|
            pileupObj = PileupLine.new(out)
            self[pileupObj.pos] = pileupObj
          end

        end
end
