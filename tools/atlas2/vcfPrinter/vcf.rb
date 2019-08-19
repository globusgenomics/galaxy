require 'vcf_variant.rb'

class Vcf < Hash

  attr_accessor :lines

  def initialize(vcfFile, *args)
    if args.length != 0
      groups = %r/(\w+):(\d+)-(\d+)/.match(args[0])
      chr = groups[1]
      start_pos = groups[2]
      end_pos = groups[3]
    end

    if File.fnmatch("*.gz", vcfFile)
      @lines = self.get_region(vcfFile, chr, start_pos, end_pos).strip.split("\n")
    else
      @lines = File.open(vcfFile,'r').readlines
    end

    @lines.each do |line|
      next if line =~ /^#/
      vcfObj = VcfLine.new(line)
      self[vcfObj.uniqPos] = vcfObj
    end
  end

  def get(*keys)
    keys.each do |key|
      return self[key] if self.has_key?(key)
    end
    raise KeyError
  end

  def get_region(vcfFile, chr, start_pos, end_pos)
    return `tabix #{vcfFile} #{chr}:#{start_pos.to_i}-#{end_pos.to_i}`
  end

end

