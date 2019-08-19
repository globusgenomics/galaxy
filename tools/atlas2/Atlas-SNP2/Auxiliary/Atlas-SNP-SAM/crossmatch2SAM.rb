=begin

Changelog:
02/15/2011 -Jin
fix header
remove incompatible MD tag
09/29/2010 -Jin
add compatibility to ruby 1.9


  *Name: crossmatch2SAM.rb
  *Description:
      'Crossmatch2SAM.rb' is designed to convert crossmatch output into SAM format.
  *Author: Yong (Tony) Wang, Zhengzheng Wan and Fuli Yu @ BCM - Human Genome Sequencing Center
  *Date:  07/29/2009

=end



#!/usr/bin/env ruby
require 'getoptlong'
require 'bigdecimal'
require 'zlib'

class  Crossmatch2SAM

#Initialize options 
def Crossmatch2SAM.processArguments()

opts = GetoptLong.new(
    ["--input", "-i", GetoptLong::REQUIRED_ARGUMENT],
    ["--help","-h", GetoptLong::NO_ARGUMENT],
    ["--output","-o", GetoptLong::REQUIRED_ARGUMENT],
    ["--read","-r", GetoptLong::REQUIRED_ARGUMENT],
    ["--qual","-q", GetoptLong::REQUIRED_ARGUMENT],
    ["--refer","-f", GetoptLong::REQUIRED_ARGUMENT]
)

optHash = {}
opts.each do |opt, arg|
  optHash[opt] = arg
end

Crossmatch2SAM.usage() if (optHash.key?("--help"));

Crossmatch2SAM.usage() if (optHash.empty?);
return optHash
end

#The usage information of the program
def Crossmatch2SAM.usage(msg='')
    unless (msg.empty?)
        puts "\n#{msg}\n"
    end
    puts "
    
Program Overview:

'Crossmatch2SAM_v1.rb' is designed to convert crossmatch output into SAM format.
   
Command Line Arguments:
  Usage: Crossmatch2SAM_v1.rb -i [crossmatch file] -o [SAM format outputfile] -r [read sequence file] -q [read quality file] -f [reference sequence file]
  -i,--input     crossmatch format file
  -o,--output    user-defined output file name
  -r,--read      gzip compressed fasta format query sequence file.
  -q,--qual      gzip compressed fasta format quality file.
  -f,--refer     fasta format reference genome file.
  
Example:

  ruby Crossmatch2SAM_v1.rb -i crossmatch.txt -r query.txt.gz -q quality.txt.gz -f reference.txt -o sample.sam 

Output format:
  query_name[tab]FLAG[tab]reference_name[tab]leftmost_postion[tab]mapping_quality[tab]CIGAR[tab]MRNM[tab]MPOS[tab] \
  ISIZE[tab]query_sequence[tab]read_quality[tab]program[tab]Alignment score[tab]number of perfect hits[tab] \
  string for mismatch positions[tab]sub_percentage[tab]del_percentage[tab]ins_percentage[tab]
  refer to SAM Format Specifications for detail"
  exit(2);
end

#Initialize the object of the class "Crossmatch2SAM"
def initialize(optHash)
    @optHash=optHash
    setParameters()
end

#Initialize the parameters
def setParameters()
    if (!@optHash.key?('--input') or !@optHash.key?("--output") or !@optHash.key?("--read") or !@optHash.key?("--qual") or !@optHash.key?("--refer"))then
        Crossmatch2SAM.usage("Option missing!")
        exit(2);
    end
    
    @crossFile=@optHash["--input"]
    @outputFile=@optHash["--output"]
    @queryFile=@optHash["--read"]
    @qualityFile=@optHash["--qual"]
    @referenceFile=@optHash["--refer"]
        
end

$pattern=/^\s*(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\S+)\s+([0-9]+)\s+(\d+)\s+\((\d+)\)\s+(C*?)\s*?(\S+)_(\d+)_(\d+)\s+\(*?(\d+)\)*?\s+(\d+)\s+\(*?(\d+)\)*?/

#Read crossmatch file, query sequence file, query quality file and reference file. Then write the information to output file in SAM format
def readcrossmatchwritesam()
    
    #control parameters
    $start=0
    $found=0
    $num=0
    $num_id=0
    $pre_refer = ""
    $cur_refer = ""
    $prev_num = 0
    $prev_type = ""

    #output information parameters
    #Alignment section parameters
    $qseq = ""
    $qualseq = ""
    $cigar=""
    type = []
    type_id = []
    type_num = []
    type_num_id = []
    qplace = []
    qplace_id = []
    tplace = []
    tplace_id = []
    ntype = []
    nqual = []
    nseq = []
    qualityArray = []
    referArray = []
    md_tag_array = []

    $lengthQuery = 0
    $lengthCross = 0
    $trim = 0

    #Tag section parameters
    $tag_NM =0
    $tag_H0 = 0
#    $tag_MD = ""
    $tag_C = 0

    #Get file handlers
    crossReader=File.open(@crossFile,"r")
    samWriter=File.new(@outputFile,"w")
    #queryReader=File.open(@queryFile,"r")
    queryReader=Zlib::GzipReader.open(@queryFile)
    #qualityReader=File.open(@qualityFile,"r")
    qualityReader=Zlib::GzipReader.open(@qualityFile)
    referenceReader=File.open(@referenceFile,"r")
   
    #read query sequence file into a read hash table
    qseqHash = Hash.new
    $queryname=""
    queryReader.each do |qline|
        if qline =~ /^>/
            $queryname = qline.split(/[ \t>\n]+/)[1] 
            qseqHash[$queryname] = ""
        else
            qseqHash[$queryname] << qline.chomp.strip
        end
    end 
    $stderr.puts "Good: finished hashing read file"

    #read quality file into a quality hash table
    qualityHash = Hash.new
    #open the quality sequence file and then convert it to Phred base score format (MAQ)  
    $qualseq = ""
    $queryname = ""
    qualityReader.each do |qline|
        if qline =~ /^>/
            $queryname = qline.split(/[ \t>\n]+/)[1] 
            qualityHash[$queryname] = ""
        else
            $qualseq = qline.chomp.strip
            qualityArray = $qualseq.split(/\s* \s*/)
            $qualseq = ""
            qualityArray.length.times do |i|
                temp=qualityArray[i].to_i+33
                $qualseq << temp.chr
            end
            qualityHash[$queryname] << $qualseq
        end
    end
    $stderr.puts "Good: finished hashing quality file"  

    #read reference sequence file into hash table
    $referenceHash = Hash.new
    $referenceName=""
    referenceReader.each do |qline|
        if qline =~ /^>(\S+)/
            $referenceName = $1
	    if $referenceName =~ /^chr(\S+)/ 
                $referenceName = $1
            end
            $referenceHash[$referenceName]=""
        else
            $referenceHash[$referenceName] << qline.chomp.strip
        end
    end
    #puts $referenceHash[$referenceName]
    $stderr.puts "Good: finished hashing reference file"
#    $sequenceLength = $referenceHash[$referenceName].length

    #print header information
    samWriter.print "@HD\tVN:1.0\n"
    $referenceHash.keys.each do |ref|
        samWriter.print "@SQ\tSN:#{ref}\tLN:#{$referenceHash[ref].length}\n"   
    end
#    samWriter.print "@SQ\tSN:#{$referenceName}\tLN:#{$sequenceLength}\n"
    samWriter.print "@RG\tID:--\tSM:--\tCN:HGSC\n"
    samWriter.print "@PG\tID:crossmatch\n"
                            
    crossReader.rewind
    crossReader.each do |line|

        if line.match($pattern) or crossReader.eof?
            if( crossReader.eof? )
                line=~ /^\s(\S+)\s+(\d+)\s+(\S+)\((\d+)\)\s+(\d+)\s+(\S+)/
                temp_type, qplace[$num], ntype[$num], nqual[$num], tplace[$num], nseq[$num] = $1, $2.to_i, $3, $4.to_i, $5.to_i, $6
                if temp_type =~ /^(\S+)\-(\d+)/
                    type[$num] = $1
                    type_num[$num] = $2.to_i
                else
                    type[$num] = temp_type
                    type_num[$num] = 1
                end
                $num = $num + 1
            end

            if $start!=0
             
                #get quality ASCII sequence from quality sequence hash
                $qualseq = ""
                if qualityHash.has_key?($qname)
                    $qualseq = qualityHash[$qname]
                else
                    $stderr.puts "Warning:  can't find the quality sequence for the query:\t#{$qname}\n"
                end

                #get query sequence from query sequence hash
                $qseq = ""
                if qseqHash.has_key?($qname)
                    $qseq = qseqHash[$qname]
                    #$stderr.puts $qseq
                    $lengthQuery = $qseq.length
                    $lengthQual = $qualseq.length
                    $lengthCross = $eposq + $nofbases
                    temp = $lengthQuery - $lengthCross
                    
                    if temp >0 
                        #read length larger than cross length, need to trim the read sequence
                        if temp <= 20
                             #trim the head adapter
                             $qseq = $qseq[temp..-1]
                        else
                             #Trim the head adapter and tail adapter
                             $qseq = $qseq[20..-1]
                             temp = 19 - temp
                             $qseq = $qseq[0..temp]
                        end
                    elsif temp == 0
                        #perfect match, nothing needed to do
                    else
                        #read length smaller than cross length, potential problem with read sequence
                        $stderr.puts "read length smaller than cross length\n"
                        #exit(0);
                    end
 
                    temp = $lengthQual - $lengthCross
                    if temp >0
                        #read qual length larger than cross length, need to trim the read qual sequence
                        if temp <= 20
                             #trim the head adapter
                             $qualseq = $qualseq[temp..-1]
                        else
                             #Trim the head adapter and tail adapter
                             $qualseq = $qualseq[20..-1]
                             temp = 19 - temp
                             $qualseq = $qualseq[0..temp]
                        end
                    elsif temp == 0
                        #perfect match, nothing needed to do
                    else
                        #read length smaller than cross length, potential problem with read sequence
                        $stderr.puts "read length smaller than cross length\n"
                        #exit(0);
                    end

                else
                    $stderr.puts "Warning:  can't find the query sequence for the query:\t#{$qname}\n"
                end

                $tag_NM = 0
                if $compl == "C"
                    if $sposr < $eposr
                       $pos = $pos1 + $sposr
                       $flag=0
                    else
                       $pos = $pos1 + $eposr
                       $flag=16
                    end
                else
                    $flag=0
                    if $sposr < $unk
                       $pos = $pos1 + $sposr
                    else
                       $pos = $pos1 + $unk
                    end
                end
                $pos -= 1

                if $compl == "C"
                    $tag_C = 0
                else
                    $tag_C = 1
                end

                               
                if $num > 0

                    #Start to formulate CIGAR code

                    #Count the number of positions of I and D, since S is counted as mismatch
                    #Also count the number of nucleotide differences and save to tag_NM
                    $num_id = 0
                    (0..$num-1).each do |i|
                        $tag_NM += type_num[i]
                        if ( type[i] == "I" or type[i] == "D")
                            type_id[$num_id]=type[i]
                            type_num_id[$num_id]=type_num[i]
                            qplace_id[$num_id]=qplace[i]
                            $num_id = $num_id + 1
                        end 
                    end

                    #Use the number of positions of I and D to formulate CIGAR code
                    
                    if $num_id > 0
                        #If there are I or D
                        $prev_num = 1
                        if $sposq > 1
                            #Soft clips in the beginning
                            temp = $sposq-1 
                            $cigar = temp.to_s + "S"
                            #$stderr.puts "softclips #{$qname}\n"
                        else
                            $cigar = ""
                        end

                        if type_id[0] == "I"
                            #The first is "I"
                            temp = qplace_id[0] - $sposq 
                        else
                            #The first is "D"
                            temp = qplace_id[0] - $sposq + 1
                        end
                        if temp > 0
                            $cigar << temp.to_s + "M"
                        end
                        $cigar << type_num_id[0].to_s + type_id[0]
                        $prev_num = type_num_id[0]
                        $prev_type = type_id[0]
                        
                        #Calculate "I" and "D"
                        (1..$num_id-1).each do |i|
                            if type_id[i] == "I"
                                if $prev_type == "I"
                                    temp = qplace_id[i] - qplace_id[i-1] - $prev_num
                                else
                                    temp = qplace_id[i] - qplace_id[i-1] - 1
                                end
                            else
                                if $prev_type == "I"
                                    temp = qplace_id[i] - qplace_id[i-1] - $prev_num + 1
                                else
                                    temp = qplace_id[i] - qplace_id[i-1]
                                end  
                            end
                            $prev_num = type_num_id[i]
                            $prev_type = type_id[i]

                            if temp > 0
                                $cigar << temp.to_s + "M"
                            end
                            
                            $cigar << type_num_id[i].to_s + type_id[i]
                        end

                        #Calculate the last one
                        if type_id[$num_id-1] =~ /^(\S+)\-(\d+)/
                            if $1 == "I"
                                temp = $eposq - qplace_id[$num_id-1] - $2.to_i + 1 
                            else
                                temp = $eposq - qplace_id[$num_id-1] 
                            end
                        else
                            if type_id[$num_id-1] == "I"
                                temp = $eposq - qplace_id[$num_id-1]
                            else
                                temp = $eposq - qplace_id[$num_id-1]
                            end
                        end
                        if type_num_id[$num_id-1] > 1
                            if type_id[$num_id-1] == "I"
                                temp = $eposq - qplace_id[$num_id-1] - type_num_id[$num_id-1] + 1 
                            else
                                temp = $eposq - qplace_id[$num_id-1] 
                            end
                        else
                            temp = $eposq - qplace_id[$num_id-1]
                        end
                        if temp > 0
                            $cigar << temp.to_s + "M"
                        end
                        
                     else
                        #If there is no I or D, all are mismatch
                        if $sposq > 1
                            #Soft clips in the beginning
                            temp = $sposq-1 
                            $cigar = temp.to_s + "S"
			    temp = $eposq - $sposq + 1
                            if temp > 0
                                $cigar << temp.to_s + "M"
                            end
                            
                        else
                            temp = $eposq - $sposq + 1
                            if temp > 0
                                $cigar = temp.to_s 
                                $cigar << "M"
                            end
                        end
     
                     end
                     if $nofbases > 0
                        #Soft clips in the end
                        $cigar << $nofbases.to_s + "S"
                     end
		     #puts $cigar
                     #CIGAR code is formulated and ready to output
                     
                     if $flag == 16
                         $cigar = $cigar.split(/(\d+[MIDNSHP])/).reverse.join
                         $qseq = $qseq.reverse
                         $qseq = $qseq.tr('ACGT','TGCA')
                         $qualseq = $qualseq.reverse
                     end
#                     $tag_MD = getMDtag()

                     samWriter.print "#{$qname}\t#{$flag}\t#{$rname}\t#{$pos}\t99\t#{$cigar}\t*\t0\t0\t#{$qseq}\t#{$qualseq}\tPG:Z:crossmatch\tAS:i:#{$score}\tNM:i:#{$tag_NM}\tH0:i:0\tXP:f:#{$psub}\tXD:f:#{$pdel}\tXI:f:#{$pins}\tXS:i:#{$tag_C}\n"
                else
                     #Perfect match
                     $tag_H0 = 1
#                     $tag_MD = ""
                     $cigar = ""
                     if $sposq > 1
                         #Soft clips in the beginning
                         temp = $sposq-1
                         $cigar = temp.to_s + "S"
                         $tag_H0 = 0
                     end
                     temp = $eposq - $sposq + 1
                     if temp > 0
                         $cigar << temp.to_s + "M"
                     end
                     if $nofbases > 0
                         #Soft clips in the end
                         $tag_H0 = 0
                         $cigar << $nofbases.to_s + "S"
                     end
                     if $flag == 16 
                         $cigar = $cigar.split(/(\d+[MIDNSHP])/).reverse.join
                         $qseq = $qseq.reverse
                         $qseq = $qseq.tr('ACGT','TGCA')
                         $qualseq = $qualseq.reverse
                     end
#                     $tag_MD = getMDtag()
                     samWriter.print "#{$qname}\t#{$flag}\t#{$rname}\t#{$pos}\t99\t#{$cigar}\t*\t0\t0\t#{$qseq}\t#{$qualseq}\tPG:Z:crossmatch\tAS:i:#{$score}\tNM:i:0\tH0:i:#{$tag_H0}\tXP:f:#{$psub}\tXD:f:#{$pdel}\tXI:f:#{$pins}\tXS:i:#{$tag_C}\n"
                end
            else
                $start=1
            end
            
            $num=0
            $num_id=0
            line.match($pattern)
            $score,$psub,$pdel,$pins,$qname,$sposq,$eposq,$nofbases,$compl,$rname,$pos1, $pos2, $unk,$sposr,$eposr =$1.to_i,$2.to_f,$3.to_f,$4.to_f,$5,$6.to_i,$7.to_i,$8.to_i,$9,$10,$11.to_i, $12.to_i, $13.to_i,$14.to_i,$15.to_i
            #puts "#{$score}\t#{psub}\t#{pdel}\t#{pins}\t#{$qname}\t#{$sposq}\t#{$eposq}\t#{$nofbases}\t#{$compl}\t#{$rname}\t#{$pos1}\t#{$pos2}\t#{unk}\t#{$sposr}\t#{$eposr}\n"
            if $rname =~ /^chr(\S+)/ 
                $rname = $1
            end
            $cur_refer = $rname
               
        elsif line=~ /^\s(\S+)\s+(\d+)\s+(\S+)\((\d+)\)\s+(\d+)\s+(\S+)/
            temp_type, qplace[$num], ntype[$num], nqual[$num], tplace[$num], nseq[$num] = $1, $2.to_i, $3, $4.to_i, $5.to_i, $6
            if temp_type =~ /^(\S+)\-(\d+)/
                type[$num] = $1
                type_num[$num] = $2.to_i
            else
                type[$num] = temp_type
                type_num[$num] = 1
            end
            $num = $num + 1
        end
    end
    
    crossReader.close
    queryReader.close
    samWriter.close
    qualityReader.close
    referenceReader.close
end

def getMDtag()
    insertion=Hash.new()    
    deletion=Hash.new()
    array=$cigar.split(/(\d+)/)
    size=array.length
    size=(size-1)/2
    out_num=0
    start=0
    read_size=0
    if array[2] == "S"
        out_num += array[1].to_i
        start += array[1].to_i
        read_size += array[1].to_i
    end
    (0..(size-1)).each do |i|
        if array[i*2+2] == "M"
            out_num += array[i*2+1].to_i
            read_size += array[i*2+1].to_i
            next
        end
        if array[i*2+2] == "D"
            deletion[out_num]=array[i*2+1].to_i
        end
        if array[i*2+2] == "I"
            (0..(array[i*2+1].to_i - 1)).each do |j|
                insertion[out_num] = 1
                out_num += 1
            end
            read_size += array[i*2+1].to_i
        end
    end
    out_num=-1
    $MD_tag = ""
    ref_position = $pos - 1
    prev_type = "ND"
    (start..(read_size-1)).each do |i|
        out_num += 1
        next if insertion.has_key?(i)
        if deletion.has_key?(i)
            if out_num > 0
                $MD_tag << out_num.to_s
            end
            $MD_tag << "^"
            prev_type = "D"
            (0..(deletion[i] - 1)).each do |j|
                ref_position += 1
                $MD_tag << $referenceHash[$cur_refer][ref_position, 1] 
            end
            out_num=0
        end
        if $qseq[i,1].upcase != $referenceHash[$cur_refer][ref_position, 1].upcase 
            if out_num == 0
                if prev_type == "D"
                    $MD_tag << out_num.to_s
                end
            else
                $MD_tag << out_num.to_s
            end
            out_num = -1
            $MD_tag << $referenceHash[$cur_refer][ref_position, 1]
        end
        prev_type = "ND"
        ref_position += 1
    end 
    out_num += 1
    #$stderr.puts "#{$pos}\n" 
    #$stderr.puts "#{$cigar}\n" 
    #temp_seq = $referenceHash[$cur_refer][($pos-1)..($pos - 1 + read_size-1)]
    #$stderr.puts "#{temp_seq}\n"
    #temp_seq = $qseq[start..(read_size-1)]
    #$stderr.puts "#{$qseq}\n"
    #if out_num > 0
    #    $MD_tag << out_num.to_s
    #end
    return $MD_tag
end

end

#Main methods
optHash=Crossmatch2SAM.processArguments()
Crossmatch2SAMDeterminter=Crossmatch2SAM.new(optHash)
Crossmatch2SAMDeterminter.readcrossmatchwritesam()
exit(0);
