#!/usr/bin/env python
"""A script which runs jobs on the condor cluster. Requires a shared file
system for both data files and the executable.
"""
import sys
import subprocess
import time
import re
import tempfile
import os
import math

SOFTWARE_PATH="/nfs/software/blast/bin"
SCRATCH_PATH="/nfs/scratch/condor_run"

if not os.path.isdir(SCRATCH_PATH):
    os.mkdir(SCRATCH_PATH)

JOB_STATUS_UNEXPANDED=0
JOB_STATUS_IDLE=1
JOB_STATUS_RUNNING=2
JOB_STATUS_REMOVED=3
JOB_STATUS_COMPLETED=4
JOB_STATUS_HELD=5
JOB_STATUS_ERROR=6

JOB_STATUS = {
    "0": ("Unexpanded", "U"),
    "1": ("Idle ", "I"),
    "2": ("Running ", "R"),
    "3": ("Removed ", "X"),
    "4": ("Completed ", "C"),
    "5": ("Held ", "H"),
    "6": ("Submission_err ", "E")
}
TEMPLATE = """
universe = vanilla
executable = %(command)s
arguments = "%(arguments)s"
output = out
error = err
log = log
initialdir = %(jobdir)s.part.$(Process)
queue %(parts)s
"""

def _quote_arg(arg):
    # escape quotes
    arg = arg.replace("'", "''")
    arg = arg.replace('"', '""')
    # quote
    return "'%s'" % arg

_proc_re = re.compile(r"^\*\* Proc (\d+)\.(\d+):$", re.MULTILINE)
_split_re = re.compile(r";SPLIT\(([^)]+)\);(.+)")
_output_re = re.compile(r";OUTPUT;(.+)")

class CondorJob(object):
    UNSUBMITTED = "unsubmitted"
    SUBMITTED = "submitted"
    COMPLETE = "complete"

    def __init__(self, command, args):
        self.command = command
        self.split_arg = None
        self.output_file = None
        self.args = []
        self.parts = 1
        self._process_args(args)
        self.job_dir = tempfile.mkdtemp(dir=SCRATCH_PATH)
        self.outdata = None
        self.errdata = None
        self.logdata = None
        self.cluster_id = None
        self.exit_code = None
        self.state = CondorJob.UNSUBMITTED

        self.condor_job = self._make_condor_job()

    def _fasta_split(self):
        infile = self.args[self.split_arg]
        self.args[self.split_arg] = "split_in"

        size = os.path.getsize(infile)
        bytes_per_part = math.ceil(size / self.parts)

        part = 0
        with open(infile) as inf:
            b = bytes_per_part + 1
            outf = None
            for line in inf:
                if line.startswith(">"):
                    if b > bytes_per_part:
                        part += 1
                        part_dir = os.path.join(self.job_dir, "part%d" % part)
                        os.mkdir(part_dir)
                        outfile = os.path.join(part_dir, "split_in")
                        if outf:
                            outf.close()
                        outf = open(outfile, "w")
                        b = 0
                outf.write(line)
                b += len(line)

        # in case there are less parts then expected, because of the
        # need to split on sequence boundaries.
        self.parts = part

    def _process_args(self, args):
        for i, arg in enumerate(args):
            m = _split_re.match(arg)
            if m:
                if self.split_arg:
                    raise ValueError("Only one split arg is supported")
                self.split_arg = i
                self.parts = int(m.group(1))
                self.args.append(m.group(2))
            else:
                m = _output_re.match(arg)
                if m:
                    self.output_file = m.group(1)
                    self.args.append("split_out")
                else:
                    self.args.append(arg)

    def _join_output(self):
        # TODO: cat all parti/split_out > self.output_file
        pass

    def __del__(self):
        try:
            os.rmdir(self.job_dir)
        except:
            pass

    def _make_condor_job(self):
        if self.parts > 1:
            self._fasta_split()
        else:
            os.mkdir(os.path.join(self.job_dir, "part1"))
        arguments = " ".join(map(_quote_arg, self.args))
        return TEMPLATE % dict(SOFTWARE_PATH=SOFTWARE_PATH,
                               command=self.command,
                               arguments=arguments,
                               jobdir=self.job_dir,
                               parts=self.parts)

    def submit(self):
        if self.state != CondorJob.UNSUBMITTED:
            raise ValueError("This class is designed for single submits only")
        child = subprocess.Popen(["condor_submit", "-"], stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = child.communicate(input=self.condor_job)
        return_code = child.returncode
        if return_code != 0:
            raise ValueError("condor_submit failed with code %d: %s"
                             % (return_code, stderr))
        # Look for first output line '** Proc 25.X:'
        # Assume all such lines have the same major id (before the dot).
        m = _proc_re.search(stdout)
        if m is None:
            raise ValueError("No Proc line in submit output")

        self.cluster_id = m.group(1)
        self.state = CondorJob.SUBMITTED
        return self.cluster_id


    def _exit_code(self):
        if self.state != CondorJob.COMPLETE:
            raise ValueError("Can't get exit code of incomplete job")
        child = subprocess.Popen(["condor_history", self.cluster_id]
                                 + "-format %s ExitCode".split(),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        stdout, stderr = child.communicate()
        try:
            return int(stdout)
        except Exception:
            return None

    def _read_tmp(self):
        if self.state != CondorJob.COMPLETE:
            raise ValueError("Can't get output of incomplete job")
        for i in xrange(1, self.parts + 1):
            part_dir = os.path.join(self.job_dir, "part" + i)
            with open(os.path.join(part_dir, "out")):
                self.outdata += f.read()
            with open(os.path.join(part_dir, "err")):
                self.errdata += f.read()
            with open(os.path.join(part_dir, "log")):
                self.logdata += f.read()

    def wait(self, sleep_time=10):
        if self.state == CondorJob.COMPLETE:
            return self.exit_code
        elif self.state != CondorJob.SUBMITTED:
            raise ValueError("Can't wait on unsubmitted job")

        while True:
            child = subprocess.Popen(["condor_q", self.cluster_id]
                                     + "-format %s JobStatus".split(),
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
            stdout, stderr = child.communicate()
            # When the job is complete, condor_q will output nothing - the
            # details are then available in condor_history. Not sure what
            # happens when the job is held; for now just assume it's
            # complete when the output is empty.

            # TODO: when it's held, it stays in the q with status 5;
            # should handle this case. Happened when my local cluster
            # was not setup correctly - the job got automatically held.
            if len(stdout) == 0:
                break
            time.sleep(sleep_time)
            del child
        self.state = CondorJob.COMPLETE
        self.exit_code = self._exit_code()
        self._read_tmp()
        self._join_output()
        # TODO: concatenate output
        return self.exit_code

    def __str__(self):
        return self.condor_job


def main():
    if len(sys.argv) < 2:
        raise ValueError("At least one argument is required")
    command = sys.argv[1]
    args = sys.argv[2:]
    condor_job = CondorJob(command, args)
    cluster_id = condor_job.submit()
    print "Job submitted with cluster_id %s" % cluster_id
    result = condor_job.wait()
    if condor_job.logdata:
        sys.stdout.write(condor_job.logdata)
    if condor_job.outdata:
        sys.stdout.write(condor_job.outdata)
    if condor_job.errdata:
        # galaxy looks for output to stderr to determine "failure",
        # not the exit code, so we send stderr to stdout unless we got
        # a failure exit code from blast.
        if result == 0:
            sys.stdout.write(condor_job.errdata)
        else:
            sys.stderr.write(condor_job.errdata)

result = 1

try:
    cmd = " ".join(sys.argv[1:])
    result = main()
except Exception, err:
    import traceback
    sys.stderr.write("Error invoking command:\n%s\n" % cmd)
    traceback.print_exc()

sys.exit(result)
