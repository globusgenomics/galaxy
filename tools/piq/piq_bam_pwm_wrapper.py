#!/usr/bin/python

import optparse, os, re, shutil, sys, tempfile, glob, shlex, vcf, pysam, tarfile, csv, operator
from subprocess import *
import subprocess
import multiprocessing
from joblib import Parallel, delayed
from itertools import izip
from os.path import isfile, join; from os import listdir

NCPU = multiprocessing.cpu_count()
CHUNK_SIZE = 2**20 #1mb

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

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

def parallel_jobs(i, piq_dir_path, motif_dir, tmp_dir, out_dir, output_dir, bamRData):
    motif=int(i)
    cmd = "Rscript %s/pertf.r %s/common.r %s/ %s %s/ %s %d; "  % (piq_dir_path, piq_dir_path, motif_dir, tmp_dir, out_dir, bamRData, motif)
#    piq_cmd = piq_cmd2+piq_cmd3
#    cmd.append(piq_cmd)
#    cmd="".join(cmd)
    print "cmd:%s" % cmd
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
    parser.add_option( '', '--motif-dir', dest='motif_dir', action='store', type="string", default=None, help='If specified, the output directory for extra files.' )

    (options, args) = parser.parse_args()

    ephemeral_dir = tempfile.mkdtemp(prefix="piq-")
    ephemeral_out_dir = "%s/out" % ephemeral_dir;
    if not os.path.exists(ephemeral_out_dir):
        os.mkdir(ephemeral_out_dir)
    ephemeral_tmp_dir = "%s/tmp" % ephemeral_dir;
    if not os.path.exists(ephemeral_tmp_dir):
        os.mkdir(ephemeral_tmp_dir)

    # get the max number of motif
    onlyfiles = [f for f in listdir(options.motif_dir) if isfile(join(options.motif_dir, f))]

    sortedfiles=[x for x in sorted_nicely(onlyfiles)]
    maxMotif=int(sortedfiles[-1].split('.')[0])

    Parallel(n_jobs=NCPU)(delayed(parallel_jobs)(i, piq_dir_path, options.motif_dir, ephemeral_tmp_dir, ephemeral_out_dir, ephemeral_dir, options.inputF) for i in range(1,maxMotif+1))
    
    #combine files
    csvFiles =sorted(glob.glob("%s/*calls.all.csv" % ephemeral_out_dir))
    csvList=create_dict(csvFiles) 
    
    for v in csvList.values():
        if len(v) == 2:
            file_merge(v[0],v[1],ephemeral_out_dir)
        else:
            outF=os.path.splitext(os.path.basename(v[0]))[0]
            outExt=os.path.splitext(os.path.basename(v[0]))[1]
            shutil.copy(v[0],"%s/%s_merged%s" % (ephemeral_out_dir, outF, outExt))

    bedFiles=sorted(glob.glob("%s/*calls.all.bed" % ephemeral_out_dir))
    bedList=create_dict(bedFiles)

    for v in bedList.values():
        if len(v) == 2:
            file_merge(v[0],v[1],ephemeral_out_dir)
        else:
            outF=os.path.splitext(os.path.basename(v[0]))[0]
            outExt=os.path.splitext(os.path.basename(v[0]))[1]
            shutil.copy(v[0],"%s/%s_merged%s" % (ephemeral_out_dir, outF, outExt))

    m_csvFiles=sorted(glob.glob("%s/*_merged.csv" % ephemeral_out_dir))
    m_bedFiles=sorted(glob.glob("%s/*_merged.bed" % ephemeral_out_dir))

    for i in range(len(m_csvFiles)):
        outF=os.path.splitext(os.path.basename(m_csvFiles[i]))[0]
        #print "\noutF: %s\n" % outF
        with open('%s/%s_all_tmp' % (ephemeral_out_dir,outF), 'wb') as res, open(m_bedFiles[i]) as f1, open(m_csvFiles[i]) as f2:
            for line1, line2 in zip(f1, f2):
                res.write("{},{}\n".format(line1.rstrip(), line2.rstrip()))
        with open('%s/%s_all_tmp' % (ephemeral_out_dir,outF), "rb") as f3:
            rdr=csv.reader(f3)

            with open('%s/%s_final.bed' % (ephemeral_out_dir,outF),'wb') as result:
                wtr = csv.writer(result, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
                for r in rdr:
                    if len(r) > 12:
                        wtr.writerow((r[0],r[1],r[2],r[3],r[5],r[9].replace('\\',""),r[10].replace('\\',""),r[11].replace('\\',""),r[12]))
    filenames=sorted(glob.glob("%s/*_final.bed" % ephemeral_out_dir))
    with open('%s/all_motif_comb_final.bed' % ephemeral_out_dir, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

    shutil.copy('%s/all_motif_comb_final.bed' % ephemeral_out_dir, options.outputF)

if __name__=="__main__": __main__()

