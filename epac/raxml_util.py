#!/usr/bin/env python

import os
import sys
import glob
import shutil
import datetime
from subprocess import call,STDOUT
from json_util import EpaJsonParser

class FileUtils:

    @staticmethod    
    def normalize_dir(dir_str):
        if dir_str and not dir_str.endswith("/"):
           dir_str += "/"
        return dir_str

    @staticmethod    
    def remove_if_exists(fname):
        try:
            os.remove(fname)
        except:
            pass

    @staticmethod    
    def rebase(fname, old_basedir, new_basedir):
        return fname.replace(old_basedir, new_basedir)    

class RaxmlWrapper:

    def __init__(self, config): 
        self.config = config
    
    def make_raxml_fname(self, stem, job_name, absolute=True):
        fname = "RAxML_" + stem + "." + job_name
        if absolute:
            return self.config.raxml_outdir + fname
        else:
            return fname            

    def make_raxml_wildcard(self, job_name):
        return self.make_raxml_fname("*", job_name)

    def cleanup(self, job_name):
        raxml_out_mask = self.make_raxml_wildcard(job_name)
        for fl in glob.glob(raxml_out_mask):
            os.remove(fl)

    def reduce_alignment(self, align_fname, job_name="reduce"):
        reduced_fname = align_fname + ".reduced"
        FileUtils.remove_if_exists(reduced_fname)
        self.run(job_name, ["-f", "c", "-s", align_fname])
        self.cleanup(job_name)
        if os.path.isfile(reduced_fname):
            return reduced_fname
        else:
            return align_fname

    def run_epa(self, job_name, align_fname, reftree_fname, optmod_fname="", silent=True):
        raxml_params = ["-f", "v", "-s", align_fname, "-t", reftree_fname]        
        if self.config.epa_use_heuristic:
            raxml_params += ["-G", str(self.config.epa_heur_rate)]
        if self.config.epa_load_optmod and optmod_fname:
            if os.path.isfile(optmod_fname):
                raxml_params += ["-R", optmod_fname]
            else:
                print "WARNING: Binary model file not found: %s" % optmod_fname
                print "WARNING: Model parameters will be estimated by RAxML"
                
        self.run(job_name, raxml_params, silent)
        
        jp_fname = self.make_raxml_fname("portableTree", job_name) + ".jplace"
        if os.path.isfile(jp_fname):
            jp = EpaJsonParser(jp_fname)
            return jp
        else:
            print "RAxML EPA run failed, please examine the log for details:\n %s" \
                    % self.make_raxml_fname("output", job_name)
            sys.exit()

    def run(self, job_name, params, silent=True):
        self.cleanup(job_name)        
        params += ["-m", self.config.raxml_model, "-n", job_name]

        if self.config.run_on_cluster:
            self.run_cluster(params)
            return;        

        if self.config.raxml_remote_call:
            call_str = ["ssh", self.config.raxml_remote_host]
        else:
            call_str = []
        call_str += self.config.raxml_cmd + params
        if silent:        
            print ' '.join(call_str) + "\n"
            out_fname = self.make_raxml_fname("output", job_name)
            with open(out_fname, "w") as fout:
                call(call_str, stdout=fout, stderr=STDOUT)
        else:        
            call(call_str)

    def run_cluster(self, params):
        if self.config.raxml_remote_call:
            qsub_call_str = ["ssh", self.config.raxml_remote_host]
        else:
            qsub_call_str = []
        
        raxml_call_cmd = self.config.raxml_cmd + params        
        for i in range(len(raxml_call_cmd)):
            if isinstance(raxml_call_cmd[i], basestring):
                raxml_call_cmd[i] = FileUtils.rebase(raxml_call_cmd[i], self.config.epatax_home, self.config.cluster_epatax_home)
        raxml_call_str = ' '.join(raxml_call_cmd)
                
        script_fname = self.config.tmp_fname("%NAME%_sub.sh")
        FileUtils.remove_if_exists(script_fname)
        shutil.copy(self.config.cluster_qsub_script, script_fname)
        qsub_job_name = "epa"        
        with open(script_fname, "a") as fout:
            fout.write("#$ -N %s\n" % qsub_job_name)
            fout.write("\n")            
            fout.write(raxml_call_str + "\n")

        cluster_script_fname = FileUtils.rebase(script_fname, self.config.epatax_home, self.config.cluster_epatax_home)
        qsub_call_str += ["qsub", "-sync", "y", cluster_script_fname]

        print raxml_call_str + "\n"
        print ' '.join(qsub_call_str) + "\n"
#        sys.exit()

        call(qsub_call_str)
        if not self.config.debug:
            FileUtils.remove_if_exists(script_fname)

    def result_exists(self, job_name):
        if os.path.isfile(self.make_raxml_fname("result", job_name)):
            return True
        else:
            return False

    def epa_result_exists(self, job_name):
        if os.path.isfile(self.make_raxml_fname("labelledTree", job_name)):
            return True
        else:
            return False

    def copy_result_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("result", job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_optmod_params(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("binaryModelParameters", job_name)
        shutil.copy(src_fname, dst_fname)

    def copy_epa_orig_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("originalLabelledTree", job_name)
        shutil.copy(src_fname, dst_fname)
        
    def copy_epa_result_tree(self, job_name, dst_fname):
        src_fname = self.make_raxml_fname("labelledTree", job_name)
        shutil.copy(src_fname, dst_fname)
