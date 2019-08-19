#!/usr/bin/ruby
require 'settings.rb'

class ClusterJobStatus
	attr_reader :jobs

	def initialize
		rawOut = self.getRunningJobs
		@jobids={'Q'=>[], 'R'=>[], 'H'=>[], 'W'=>[]}
		@jobs=Array.new
		rawOut.each do |out|
			if out =~ /^[0-9]/
				out=out.split()
				@jobids[out[9]].push(out[0].gsub('.sug-moab', '')) if @jobids.has_key? out[9] 
				@jobs.push(out[0].gsub('.sug-moab', '')) if @jobids.has_key? out[9]
			else
				next
			end
		end
	end

        def getRunningJobs
            #out = `bjobs -u $USER`.strip.split("\n")
            out = `qstat -u #{USER}`.strip.split("\n") ######EDIT######
            if out[0] =~ /Unable to query/
                self.getRunningJobs
            end
            return out
        end

	def activeJobs
            return @jobids['R']
	end

	def eligibleJobs
            return @jobids['Q']
	end

end

#Testing Class ClusterJobStatus
if __FILE__==$0
	clusterObj=ClusterJobStatus.new
	puts "=====>Active Jobs<====="
	puts clusterObj.activeJobs
	puts "=====>Pending Jobs<====="
	puts clusterObj.eligibleJobs
end
