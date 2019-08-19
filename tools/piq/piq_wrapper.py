#!/usr/bin/python

import optparse, os, re, shutil, sys, tempfile, glob, shlex, vcf, pysam, tarfile, csv, operator
from subprocess import *
import subprocess
import multiprocessing
from joblib import Parallel, delayed
from itertools import izip

NCPU = multiprocessing.cpu_count()
CHUNK_SIZE = 2**20 #1mb

def create_dict(files):
    tmp=dict()
    t1 = [os.path.basename(i) for i in files]
    t2 = ['-'.join(i.split('-')[0:2]) for i in t1 if 'RC' not in i]
    for i in t2:
        for j in files:
            if i in j and i in tmp:
                tmp[i].append(j)
            elif i in j and i not in tmp:
                tmp[i] = [j]
            else:
                continue
    return tmp

def file_merge(file1, file2, out_dir):
    listF = (file1, file2)
    outF=os.path.splitext(os.path.basename(file1))[0]
    outExt = os.path.splitext(os.path.basename(file1))[1]
    with open("%s/%s_merged%s" % (out_dir,outF,outExt),"wb") as fw:
        for i in range(len(listF)):
            with open(listF[i],"rb") as f:
                writer = csv.writer(fw, delimiter=',', quotechar='', quoting=csv.QUOTE_NONE, escapechar='\\')
                for i in range(len(listF)):
                    with open(listF[i],"rb") as f:
                        file_f=csv.reader(f, delimiter='\t')
                        file_f.next()
                        for line in file_f:
                            writer.writerow(line)


def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def parallel_jobs(i, piq_dir_path, tmp_dir, inputF, out_dir, output_dir):
    cmd=[]
    motif=i
    piq_cmd1 = "Rscript %s/pwmmatch.exact.r %s/common.r %s/new_motifs2.txt %d %s/; " % (piq_dir_path, piq_dir_path, piq_dir_path, motif, tmp_dir)
#    piq_cmd1 = "Rscript %s/pwmmatch.exact.r %s/common.r %s/jasparVertebrateCore2016-519-matrices.txt %d %s/; " % (piq_dir_path, piq_dir_path, piq_dir_path, motif, tmp_dir)
    #piq_cmd1 = "Rscript %s/pwmmatch.exact.r %s/common.r %s/jasparfix.txt %d %s/; " % (piq_dir_path, piq_dir_path, piq_dir_path, motif, tmp_dir)
    piq_cmd2 = "Rscript %s/bam2rdata.r %s/common.r %s/bam.RData %s; " % (piq_dir_path, piq_dir_path, tmp_dir, inputF)
    piq_cmd3 = "Rscript %s/pertf.r %s/common.r %s/ %s %s/ %s/bam.RData %d; "  % (piq_dir_path, piq_dir_path, tmp_dir, tmp_dir, out_dir, tmp_dir, motif)
    piq_cmd = piq_cmd1+piq_cmd2+piq_cmd3
    cmd.append(piq_cmd)
    cmd="".join(cmd)
    #print "cmd:%s" % cmd
    stdout = tempfile.NamedTemporaryFile( prefix="piq-stdout-", dir=tmp_dir )
    stderr = tempfile.NamedTemporaryFile( prefix="piq-stderr-", dir=tmp_dir )
    proc = subprocess.Popen( args=cmd, stdout=stdout, stderr=stderr, shell=True, cwd=output_dir )
    return_code = proc.wait()

    if return_code:
        stderr_target = sys.stderr
    else:
        stderr_target = sys.stdout
        stderr.flush()
        stderr.seek(0)
    while True:
        chunk = stderr.read( CHUNK_SIZE )
        if chunk:
            stderr_target.write( chunk )
        else:
            break
    stderr.close()
    stdout.close()
    
def __main__():
    piq_dir_path="/mnt/galaxyTools/tools/piq/1.3/thashim"
    parser = optparse.OptionParser()
    parser.add_option('', '--input', dest="inputF", action='store', type="string", default=None, help='')
    parser.add_option( '', '--output', dest='outputF', action='store', type="string", default=None, help='')
    ##parser.add_option( '', '--out-dir', dest='output_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    #if not os.path.exists(options.output_dir): 
    #    os.mkdir(options.output_dir)
    #out_dir = "%s/out" % options.output_dir; 
    #if not os.path.exists(out_dir):
    #    os.mkdir(out_dir)
    #tmp_dir = "%s/tmp" % options.output_dir; 
    #if not os.path.exists(tmp_dir):
    #    os.mkdir(tmp_dir)

    ephemeral_dir = tempfile.mkdtemp(prefix="piq-")
    ephemeral_out_dir = "%s/out" % ephemeral_dir;
    if not os.path.exists(ephemeral_out_dir):
        os.mkdir(ephemeral_out_dir)
    ephemeral_tmp_dir = "%s/tmp" % ephemeral_dir;
    if not os.path.exists(ephemeral_tmp_dir):
        os.mkdir(ephemeral_tmp_dir)

    maxMotif=2
#    maxMotif=519
#    Parallel(n_jobs=NCPU)(delayed(parallel_jobs)(i, piq_dir_path, tmp_dir, options.inputF, out_dir, options.output_dir) for i in range(1,maxMotif+1))
    Parallel(n_jobs=NCPU)(delayed(parallel_jobs)(i, piq_dir_path, ephemeral_tmp_dir, options.inputF, ephemeral_out_dir, ephemeral_dir) for i in range(1,maxMotif+1))
    
    #combine files
    csvFiles =sorted(glob.glob("%s/*calls.all.csv" % ephemeral_out_dir))
#    csvFiles =sorted(glob.glob("%s/*calls.all.csv" % out_dir))
    csvList=create_dict(csvFiles) 
    
    for v in csvList.values():
        if len(v) == 2:
            file_merge(v[0],v[1],ephemeral_out_dir)
#            file_merge(v[0],v[1],out_dir)
        else:
            outF=os.path.splitext(os.path.basename(v[0]))[0]
            outExt=os.path.splitext(os.path.basename(v[0]))[1]
            shutil.copy(v[0],"%s/%s_merged%s" % (ephemeral_out_dir, outF, outExt))
#            shutil.copy(v[0],"%s/%s_merged%s" % (out_dir, outF, outExt))

    bedFiles=sorted(glob.glob("%s/*calls.all.bed" % ephemeral_out_dir))
#    bedFiles=sorted(glob.glob("%s/*calls.all.bed" % out_dir))
    bedList=create_dict(bedFiles)

    for v in bedList.values():
        if len(v) == 2:
            file_merge(v[0],v[1],ephemeral_out_dir)
#            file_merge(v[0],v[1],out_dir)
        else:
            outF=os.path.splitext(os.path.basename(v[0]))[0]
            outExt=os.path.splitext(os.path.basename(v[0]))[1]
            shutil.copy(v[0],"%s/%s_merged%s" % (ephemeral_out_dir, outF, outExt))
#            shutil.copy(v[0],"%s/%s_merged%s" % (out_dir, outF, outExt))

    m_csvFiles=sorted(glob.glob("%s/*_merged.csv" % ephemeral_out_dir))
    m_bedFiles=sorted(glob.glob("%s/*_merged.bed" % ephemeral_out_dir))
#    m_csvFiles=sorted(glob.glob("%s/*_merged.csv" % out_dir))
#    m_bedFiles=sorted(glob.glob("%s/*_merged.bed" % out_dir))

    for i in range(len(m_csvFiles)):
        outF=os.path.splitext(os.path.basename(m_csvFiles[i]))[0]
        #print "\noutF: %s\n" % outF
        with open('%s/%s_all_tmp' % (ephemeral_out_dir,outF), 'wb') as res, open(m_bedFiles[i]) as f1, open(m_csvFiles[i]) as f2:
#        with open('%s/%s_all_tmp' % (out_dir,outF), 'wb') as res, open(m_bedFiles[i]) as f1, open(m_csvFiles[i]) as f2:
            for line1, line2 in zip(f1, f2):
                res.write("{},{}\n".format(line1.rstrip(), line2.rstrip()))
#        with open('%s/%s_all_tmp' % (out_dir,outF), "rb") as f3:
        with open('%s/%s_all_tmp' % (ephemeral_out_dir,outF), "rb") as f3:
            rdr=csv.reader(f3)

            with open('%s/%s_final.bed' % (ephemeral_out_dir,outF),'wb') as result:
#            with open('%s/%s_final.bed' % (out_dir,outF),'wb') as result:
                wtr = csv.writer(result, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
                #for r in sorted_list:
                for r in rdr:
                    if len(r) > 12:
                        wtr.writerow((r[0],r[1],r[2],r[3],r[5],r[9].replace('\\',""),r[10].replace('\\',""),r[11].replace('\\',""),r[12]))
    filenames=sorted(glob.glob("%s/*_final.bed" % ephemeral_out_dir))
#    filenames=sorted(glob.glob("%s/*_final.bed" % out_dir))
    #shutil.copy('%s/*_final.bed' % out_dir, options.output_dir)
    with open('%s/all_motif_comb_final.bed' % ephemeral_out_dir, 'w') as outfile:
#    with open('%s/all_motif_comb_final.bed' % options.output_dir, 'w') as outfile:
        for fname in filenames:
#            shutil.copy('%s' % fname, options.output_dir)
#            shutil.copy('%s' % fname, ephemeral_out_dir)
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    shutil.copy('%s/all_motif_comb_final.bed' % ephemeral_out_dir, options.outputF)
#    shutil.copy('%s/all_motif_comb_final.bed' % options.output_dir, options.outputF)

    # concatenerate calls.all and RC.calls.all for csv and bed files
    #cleanup_before_exit( options.output_dir )

if __name__=="__main__": __main__()

