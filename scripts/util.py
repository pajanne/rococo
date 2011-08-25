#!/usr/bin/env python
# encoding: utf-8
'''
Created on Nov 23, 2009
by
@author: Anne Pajon (ap12)
Copyright (c) 2009 Wellcome Trust Sanger Institute. All rights reserved.
'''

import os
import subprocess

### ---------------------------------------------------------------------------
### Logging configuration
### ---------------------------------------------------------------------------
import logging
log = logging.getLogger("genepy.util")


### ---------------------------------------------------------------------------
### Utility methods
### ---------------------------------------------------------------------------
def checkFile(filename):
    """
    Check if a given filename exists and if it is a file.
    """
    if not filename:
        raise UtilException("No filename given")
    if not os.path.exists(filename):
        raise UtilException("file %s does not exist" % filename)
    else:
        if not os.path.isfile(filename):
            raise UtilException("%s is not a file" % filename)
            
### ---------------------------------------------------------------------------
def checkDir(directory):
    """
    Check if a given directory exists.
    """
    if not os.path.exists(directory):
        raise UtilException("directory %s/ does not exist." % directory)
    
### ---------------------------------------------------------------------------
def checkSoft(softname):
    """
    Check if a given software is installed, exit if not installed
    """
    if not isSoftInstalled(softname):
        raise UtilException("software %s is not installed, please install." % softname)

### ---------------------------------------------------------------------------
def isSoftInstalled(softname):
    """
    Return true if a given software is installed, false otherwise
    """
    retval = subprocess.call(['which %s' % softname], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if retval == 0:
        log.info("software %s is installed." % softname)
        return True
    else:
        log.warning("software %s is not installed, please install." % softname)
        return False

### ---------------------------------------------------------------------------
def createDir(directory):
    """
    Check if a given directory exists, and create if it does not.
    """
    #directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        raise UtilException("directory %s/ already exists, please rename or delete." % directory)

### ---------------------------------------------------------------------------
def rmDir(directory):
    """
    Remove a directory using shell command rm
    """
    if os.path.exists(directory):
        cmd = "rm -rf %s/" % directory
        runProcess(cmd)
    else:
        log.warning("directory %s/ does not exist." % directory)
        
### ---------------------------------------------------------------------------
def rmFile(file):
    """
    Remove a file 
    """
    if os.path.exists(file):
        os.remove(file)
    else:
        log.warning("file %s does not exist." % file)
        
### ---------------------------------------------------------------------------
def isLsf():
    """
    Return True if running on LSF by testing if bsub command is installed
    """
    retval = subprocess.call(['which bsub'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if retval == 0:
        log.info("Running on LSF")
        return True
    else:
        return False
    
### ---------------------------------------------------------------------------
def runProcess(cmd):
    """
    Run a command using subprocess
    """
    log.info("$ %s" % cmd)
    process = subprocess.Popen([cmd], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # communicate() if called after poll() seams to set retval to None
    process_out, process_err = process.communicate()
    retval = process.poll()
    if not retval == 0:
        raise UtilException("Error executing %s command: %s" % (cmd.split()[0], process_err))
        return None
    else:
        log.info(process_err)
        log.info(process_out)
        return process_out
    
### ---------------------------------------------------------------------------
def submitJobArray(jobname, jobnum, jobdir, cmd, maxjob=10, queue='small'):
    """
    Submit a bsub job array on LSF using subprocess
    """
    bsub_cmd = """bsub -q %s -J"%s[1-%s]%%%s" -o %s/%%J.%%I.out -e %s/%%J.%%I.err '%s'""" % (queue, jobname, jobnum, maxjob, jobdir, jobdir, cmd)
    runProcess(bsub_cmd)

### ---------------------------------------------------------------------------
def submitJob(jobname, cmd, queue='small', outdir=None):
    """
    Submit a bsub job on LSF using subprocess
    """
    if outdir:
        bsub_cmd = """bsub -q %s -J %s -o %s/_%s.o -e %s/_%s.e '%s'""" % (queue, jobname, outdir, jobname, outdir, jobname, cmd)
    else:
        bsub_cmd = """bsub -q %s -J %s -o _%s.o -e _%s.e '%s'""" % (queue, jobname, jobname, jobname, cmd)
    runProcess(bsub_cmd)

### ---------------------------------------------------------------------------
def submitJobDependency(prev_jobname, jobname, cmd, queue='small'):
    """
    Submit a bsub dependency job

    bsub -w "done(myarrayA[*])" -J "myArrayB[1-10]" myJob2
    
    bsub  -w 'dependency_expression' some_other_task
    Use the following conditions to form the dependency expression.
        done(job_ID |"job_name" ...)
            The job state is DONE.
            LSF refers to the oldest job of job_name in memory.
        ended(job_ID | "job_name")
            The job state is EXIT or DONE.
        exit(job_ID | "job_name" [,[operator] exit_code])
            The job state is EXIT, and the job's exit code satisfies the comparison test.
            If you specify an exit code with no operator, the test is for equality (== is assumed).
            If you specify only the job, any exit code satisfies the test.
    """
    bsub_cmd = """bsub -q %s -w 'done(%s)' -J %s -o %s.bsub_out -e %s.bsub_err '%s'""" % (queue, prev_jobname, jobname, jobname, jobname, cmd)
    runProcess(bsub_cmd)

### ---------------------------------------------------------------------------
def submitJobWait(jobname, queue='small'):
    """
    Submit a bsub interactive dependency job

    bsub  -w 'dependency_expression' some_other_task
    Use the following conditions to form the dependency expression.
        done(job_ID |"job_name" ...)
            The job state is DONE.
            LSF refers to the oldest job of job_name in memory.
        ended(job_ID | "job_name")
            The job state is EXIT or DONE.
        exit(job_ID | "job_name" [,[operator] exit_code])
            The job state is EXIT, and the job's exit code satisfies the comparison test.
            If you specify an exit code with no operator, the test is for equality (== is assumed).
            If you specify only the job, any exit code satisfies the test.
    """
    log.info("Waiting for job to end...")
    bsub_cmd = """bsub -q %s -I -w 'ended(%s)' 'echo job %s ended'""" % (queue, jobname, jobname)
    runProcess(bsub_cmd)

### ---------------------------------------------------------------------------
def printHelp(options):
    """
    Return True if no argument given from OptionParser (options, args) = parser.parse_args()
    """
    for value in options.__dict__.values():
        if value != None:
            return False
    return True

### ---------------------------------------------------------------------------
def splitLine(line):
    """
    Returns the two first elements of a line splitted on ||
    """
    # ! common_name||sequence_file
    line = line.strip()
    values = line.split('||')
    return values[0], values[1]

### ---------------------------------------------------------------------------
### Class
### ---------------------------------------------------------------------------
class UtilException(Exception):
    """
    Ideally, all errors should be raised with this (or a subclass of).
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
