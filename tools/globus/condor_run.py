#!/usr/bin/env python
"""
Classes for constructing, submitting, and monitoring a Condor job. It started
as a simple script to take a normal local execution and turn it into a condor
execution, and has evolved features for executing multiple instances of
a single command with different arguments.

It's designed to work under the following assuption:

* Any file arguments, for input or output, references in the command are
  on a shared filesystem.
* There is scratch space available, also on the shared filesystem.

A temp directory under the scratch path is created for each job, and each
process within the job gets it's own subdirectory for output, error, and log
files.

When executed as a main script, the first argument is taken as the command and
the rest as arguments, and the command is run on the condor cluster. stdout
and stderr from the command is echoed.

Note: the job dir is not being deleted on exit, this is a bug. It's useful to
leave the dir around for debugging, so it should be an option, but the default
should be to clean up the temp dir.
"""
import sys
import subprocess
import time
import re
import tempfile
import os

# Temporary location for stdout, stderr, and log files.
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
error = err
log = log

initialdir = %(job_dir)s/$(Process)
"""

def _quote_arg(arg):
    # escape quotes
    arg = arg.replace("'", "''")
    arg = arg.replace('"', '""')
    # quote
    return "'%s'" % arg

_proc_re = re.compile(r"^\*\* Proc (\d+)\.(\d+):$", re.MULTILINE)

class CondorQueueItem(object):
    """
    By default, output is assumed to be debug info and is agregated with
    the output from other processes within the job. If capture to a specific
    file without agregation in desired, it can be specified here.
    """
    def __init__(self, args, output=None):
        self.args = args
        self.output = output

    def __str__(self):
        arg_string = " ".join(map(_quote_arg, self.args))
        s = '\narguments = "%s"\n' % arg_string
        if self.output:
            s += 'output = %s\n' % self.output
        else:
            s += 'output = out\n'
        s += "queue\n"
        return s

class CondorJob(object):
    UNSUBMITTED = "unsubmitted"
    SUBMITTED = "submitted"
    COMPLETE = "complete"

    def __init__(self, command, arg_lists=None, items=None):
        """
        @param arg_lists: iterable that contains (or generates) one or more
                          argument lists. Each arg list gets queued as a
                          separate instance of the command.
        """
        self.command = command

        if items:
            self.items = items
        else:
            self.items = []
        if arg_lists:
            for arg_list in arg_lists:
                self.items.append(CondorQueueItem(arg_list))

        self.parts = 0
        self.job_dir = tempfile.mkdtemp(dir=SCRATCH_PATH)
        self.outdata = ""
        self.errdata = ""
        self.logdata = ""
        self.cluster_id = None
        self.exit_code = None
        self.state = CondorJob.UNSUBMITTED

        self.condor_job = self._make_condor_job()

    def __del__(self):
        try:
            os.rmdir(self.job_dir)
        except:
            pass

    def _make_condor_job(self):
        job = TEMPLATE % dict(command=self.command,
                              job_dir=self.job_dir)
        item_strings = []
        for item in self.items:
            part_dir = os.path.join(self.job_dir, str(self.parts))
            os.mkdir(part_dir)
            item_strings.append(str(item))
            self.parts += 1
        job += "\n".join(item_strings)
        return job

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
        open(os.path.join(self.job_dir, "cid_%s" % self.cluster_id),
             "w").close()
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
            stdout = stdout.strip()
            # TODO: make sure is all ints
            print "Return codes:", stdout
            for c in stdout:
                if c != "0":
                    return 1
            return 0
        except Exception:
            return None

    def _read_tmp(self):
        if self.state != CondorJob.COMPLETE:
            raise ValueError("Can't get output of incomplete job")
        for i in xrange(self.parts):
            part_dir = os.path.join(self.job_dir, str(i))
            outfile = os.path.join(part_dir, "out")
            if os.path.exists(outfile):
                with open(outfile) as f:
                    self.outdata += f.read()
            with open(os.path.join(part_dir, "err")) as f:
                self.errdata += f.read()
            with open(os.path.join(part_dir, "log")) as f:
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
        return self.exit_code

    def __str__(self):
        return self.condor_job


def main():
    if len(sys.argv) < 2:
        raise ValueError("At least one argument is required")
    command = sys.argv[1]
    args = sys.argv[2:]
    condor_job = CondorJob(command, [args])
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

if __name__ == '__main__':
    result = 1

    try:
        cmd = " ".join(sys.argv[1:])
        result = main()
    except Exception, err:
        import traceback
        sys.stderr.write("Error invoking command:\n%s\n" % cmd)
        traceback.print_exc()

    sys.exit(result)
