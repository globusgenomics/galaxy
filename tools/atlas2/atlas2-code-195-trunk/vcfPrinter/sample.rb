require 'pileup.rb'
require 'vcf.rb'

class Sample
        
        attr_reader :name, :pileup, :vcf

        def initialize(name, vcfFile, *opt_args)
          @name=name
          if opt_args.length > 0
            @vcf = Vcf.new(vcfFile, opt_args[0])
          else
            @vcf = Vcf.new(vcfFile)
          end
          @pileupExist=false
        end

        def createPileup(pileupFile)
          @pileup = Pileup.new(pileupFile, )
        end

        def hasPileup
          return @pileupExist
        end

end

