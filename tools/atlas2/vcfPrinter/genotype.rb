#!/usr/bin/ruby

class Genotype

        def initialize(format,genotypeField)
                #fields of interest: GT,VR,RR,TD,GQ
                @fields=Hash.new()
                format.split(':').each_with_index do |f,index|
                        @fields[f]=index
                end
                begin
                        @genotypeInfo=genotypeField.split(':')
                rescue
                        puts "#{@genotypeInfo}"
                        Process.exit(0)
                end
        end

        def get(fieldOfInterest)
                if @fields.include? fieldOfInterest
                        return @genotypeInfo[@fields[fieldOfInterest]]
                else
                        puts "#{fieldOfInterest} does not exist."
                        return nil
                end
        end

        def to_s
                return @genotypeInfo.join(sep=':')
        end

end
