#!/usr/bin/env python
import os
import sys
import socket
import subprocess
import csv
import re
import multiprocessing as mp
from copy import deepcopy
from argparse import ArgumentParser
from FTAnalysisJob import BatchManager
from FTNEvents import getSkimNEvents
from getVOMSProxy import getVOMSProxy


def get_nevts(skim_tag, period, short_name, dataset):
   res = dict()
   nevts = getSkimNEvents(skim_tag, period, short_name)
   if nevts>0:
      res[dataset]=nevts
   return res


def get_strargs_from_dict(dict_ent):
   argslist = []
   for xx in dict_ent:
      if dict_ent[xx] is not None:
         argslist.append("=".join([xx,str(dict_ent[xx])]))
      else:
         argslist.append(xx)
   return " ".join(argslist)


def run_single(job_args):
   BatchManager(job_args)


def run(args):
   nthreads = min(args.nthreads, mp.cpu_count())
   hostname = socket.gethostname()
   grid_proxy, grid_user = getVOMSProxy()

   allowed_sites = args.sites
   if allowed_sites is None:
      if "t2.ucsd.edu" in hostname:
         allowed_sites = "T2_US_UCSD,T2_US_Caltech,T2_US_MIT,T3_US_UCR,T3_US_Baylor,T3_US_Colorado,T3_US_NotreDame,T3_US_Cornell,T3_US_Rice,T3_US_Rutgers,T3_US_UCD,T3_US_TAMU,T3_US_TTU,T3_US_FIU,T3_US_FIT,T3_US_UMD,T3_US_OSU,T3_US_OSG,T3_US_PuertoRico"
      elif "uscms.org" in hostname:
         #allowed_sites = "T2_US_Caltech,T2_US_MIT,T3_US_UCR,T3_US_Baylor,T3_US_Colorado,T3_US_NotreDame,T3_US_Cornell,T3_US_Rice,T3_US_Rutgers,T3_US_UCD,T3_US_TAMU,T3_US_TTU,T3_US_FIU,T3_US_FIT,T3_US_UMD,T3_US_OSU,T3_US_OSG,T3_US_PuertoRico,T2_CH_CERN,T2_BE_IIHE,T2_CN_Beijing,T2_RU_IHEP,T2_BE_UCL,T2_AT_Vienna,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CSCS,T2_DE_DESY,T2_DE_RWTH,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_KR_KNU,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_RRC_KI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_UA_KIPT,T2_UK_London_Brunel,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T3_CO_Uniandes,T3_FR_IPNL,T3_GR_IASA,T3_HU_Debrecen,T3_IT_Bologna,T3_IT_Napoli,T3_IT_Perugia,T3_IT_Trieste,T3_KR_KNU,T3_MX_Cinvestav,T3_RU_FIAN,T3_TW_NCU,T3_TW_NTU_HEP,T3_UK_London_QMUL,T3_UK_SGrid_Oxford,T3_CN_PKU"
         #allowed_sites = "T2_US_Caltech,T2_US_MIT,T3_US_UCR,T3_US_Baylor,T3_US_Colorado,T3_US_NotreDame,T3_US_Cornell,T3_US_Rice,T3_US_Rutgers,T3_US_UCD,T3_US_TAMU,T3_US_TTU,T3_US_FIU,T3_US_FIT,T3_US_UMD,T3_US_OSU,T3_US_OSG,T3_US_PuertoRico,T2_CH_CERN,T2_BE_IIHE,T2_CN_Beijing,T2_RU_IHEP,T2_BE_UCL,T2_AT_Vienna,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CSCS,T2_DE_DESY,T2_DE_RWTH,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_KR_KNU,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_RRC_KI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_UA_KIPT,T2_UK_London_Brunel,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T3_CO_Uniandes,T3_FR_IPNL,T3_GR_IASA,T3_HU_Debrecen,T3_IT_Bologna,T3_IT_Napoli,T3_IT_Perugia,T3_IT_Trieste,T3_KR_KNU,T3_MX_Cinvestav,T3_RU_FIAN,T3_TW_NCU,T3_TW_NTU_HEP,T3_UK_London_QMUL,T3_UK_SGrid_Oxford,T3_CN_PKU"
         allowed_sites = "T2_US_UCSD,T2_US_Caltech,T2_US_MIT,T3_US_UCR,T3_US_Baylor,T3_US_Colorado,T3_US_NotreDame,T3_US_Cornell,T3_US_Rice,T3_US_Rutgers,T3_US_UCD,T3_US_TAMU,T3_US_UMD,T3_US_OSU,T3_US_OSG,T3_US_PuertoRico,T2_CH_CERN,T2_BE_IIHE,T2_CN_Beijing,T2_RU_IHEP,T2_BE_UCL,T2_DE_DESY,T2_DE_RWTH,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_TR_METU"
      else:
         allowed_sites = "T2_CH_CERN,T2_BE_IIHE,T2_CN_Beijing,T2_RU_IHEP,T2_BE_UCL,T2_AT_Vienna,T2_BR_SPRACE,T2_BR_UERJ,T2_CH_CSCS,T2_DE_DESY,T2_DE_RWTH,T2_ES_CIEMAT,T2_ES_IFCA,T2_FI_HIP,T2_FR_CCIN2P3,T2_FR_GRIF_IRFU,T2_FR_GRIF_LLR,T2_FR_IPHC,T2_GR_Ioannina,T2_HU_Budapest,T2_IN_TIFR,T2_IT_Bari,T2_IT_Legnaro,T2_KR_KNU,T2_PK_NCP,T2_PL_Swierk,T2_PL_Warsaw,T2_PT_NCG_Lisbon,T2_RU_INR,T2_RU_ITEP,T2_RU_JINR,T2_RU_PNPI,T2_RU_RRC_KI,T2_RU_SINP,T2_TH_CUNSTDA,T2_TR_METU,T2_UA_KIPT,T2_UK_London_Brunel,T2_UK_London_IC,T2_UK_SGrid_Bristol,T2_UK_SGrid_RALPP,T3_CO_Uniandes,T3_FR_IPNL,T3_GR_IASA,T3_HU_Debrecen,T3_IT_Bologna,T3_IT_Napoli,T3_IT_Perugia,T3_IT_Trieste,T3_KR_KNU,T3_MX_Cinvestav,T3_RU_FIAN,T3_TW_NCU,T3_TW_NTU_HEP,T3_UK_London_QMUL,T3_UK_SGrid_Oxford,T3_CN_PKU"

   batch_script = os.path.expandvars("${CMSSW_BASE}/src/tttt/test/condor_executable_analysis.sh")
   condor_site = "t2.ucsd.edu"
   condor_outdir = "/ceph/cms/store/group/tttt/Worker/" + grid_user
   jobs_cfgdir = os.path.expandvars("${CMSSW_BASE}/src/tttt/test/output/tasks/"+args.tag)
   extratarcmd=""
   if args.extra_tar is not None:
      extratarcmd="addfile={}".format(extratarcmd)

   print("Batch processing parameters:")
   print("- Grid proxy = {}".format(grid_proxy))
   print("- Grid user = {}".format(grid_user))
   print("- Jobs directory = {}".format(jobs_cfgdir))
   print("- Batch script = {}".format(batch_script))
   print("- Batch sites = {}".format(allowed_sites))
   print("- Destination site = {}".format(condor_site))
   print("- Destination directory = {}".format(condor_outdir))
   print("- Extra uploads command = {}".format(extratarcmd))

   skim_tag = None
   output_tag = None
   extra_args = dict()
   extra_args_list = []
   csvs = []
   for rg in args.extras:
      if rg.endswith(".csv"):
         csvs.append(rg)
      else:
         extra_args_list.append(rg)
         rgfs = rg.split('=')
         if len(rgfs)==2:
            extra_args[rgfs[0]] = rgfs[1]
         else:
            extra_args[rgfs[0]] = None
         if "input_tag=" in rg:
            skim_tag = rg.replace("input_tag=","")
         elif "output_tag=" in rg:
            output_tag = rg.replace("output_tag=","")

   print("General job parameters:")
   print("- csv files = {}".format(','.join(csvs)))
   print("- Arguments to the executable = '{}'".format(' '.join(extra_args_list)))
   print("- Input ntuples tag = {}".format(skim_tag))
   print("- Output ntuples tag = {}".format(output_tag))

   csv_entries = []
   for fname in csvs:
      with open(fname) as fh:
         reader = csv.DictReader(fh)
         for row in reader:
            dataset = row["dataset"]
            if dataset.strip().startswith("#"):
               continue
            short_name = row["short_name"]
            period = row["period"]
            xsec=None
            BR=None
            isData = True
            if "AODSIM" in dataset:
               xsec=row["xsec"]
               BR=row["BR"]
               isData = False

            dict_ent=dict()
            dict_ent["dataset"]=dataset
            dict_ent["short_name"]=short_name
            dict_ent["period"]=period
            if not isData:
               dict_ent["xsec"]=xsec
               dict_ent["BR"]=BR
            csv_entries.append(dict_ent)

   nevts_map = dict()
   pool = mp.Pool(nthreads)
   [ pool.apply_async(get_nevts, args=(skim_tag, csve["period"], csve["short_name"], csve["dataset"]), callback=nevts_map.update) for csve in csv_entries ]
   pool.close()
   pool.join()

   #for dset in nevts_map:
   #   print("Number of events in '{}' = {}".format(dset, nevts_map[dset]))

   data_groups = dict()
   if args.group_eras:
      for dict_ent in csv_entries:
         dset = dict_ent["dataset"]
         if "AODSIM" in dset or dset not in nevts_map:
            continue
         sname = dict_ent["short_name"]
         period = dict_ent["period"]
         mm = re.search("(?<=Run)20[1-9]*[A-Z]",dict_ent["dataset"])
         if mm:
            pp = mm.group(0)
            if pp=="2016F":
               if "_APV" in period:
                  pp = pp+"_APV"
               else:
                  pp = pp+"_NonAPV"
            ppo = "Run{}".format(pp)
            if ppo not in data_groups:
               data_groups[ppo] = [[],[],pp] # Data sets, short names, period
            data_groups[ppo][0].append(dset)
            data_groups[ppo][1].append(sname)
         else:
            raise RuntimeError("Could not find the period string match in {}".format(dset))


   #for dd in data_groups:
      #print("Acquired data group {} for data sets {}".format(dd, ",".join(data_groups[dd][0])))

   ignored_csv_dsets = []
   all_args_list = []
   for pp in data_groups:
      dict_ent = dict()
      dsets = data_groups[pp][0]
      snames = data_groups[pp][1]

      ignored_csv_dsets.extend(dsets)

      dict_ent["dataset"] = ",".join(dsets)
      dict_ent["short_name"] = ",".join(snames)
      dict_ent["period"] = data_groups[pp][2]
      dict_ent["output_file"] = pp
      nevts_total = 0
      for dset in dsets:
         nevts_total += nevts_map[dset]
      nevts_map[dict_ent["dataset"]] = nevts_total
      all_args_list.append(dict_ent)

   for dict_ent in csv_entries:
      dset = dict_ent["dataset"]
      if dset in ignored_csv_dsets or dset not in nevts_map:
         continue
      all_args_list.append(dict_ent)

   # Add event splitting
   if args.nevts>0:
      tmp_all_args_list = []
      for dict_ent in all_args_list:
         dset = dict_ent["dataset"]
         nevts = nevts_map[dset]
         nchunks = int(nevts / args.nevts)
         if nchunks*args.nevts<nevts:
            nchunks = nchunks + 1
         if "AODSIM" not in dset:
            period = dict_ent["period"]
            runlumilist = subprocess.check_output(["print_runinfo period={}".format(period)], shell=True)
            runlumilist = runlumilist.split()
            nchunks = min(nchunks,len(runlumilist))
         for ichunk in range(0, max(1, nchunks)):
            tmp_args = deepcopy(dict_ent)
            if nchunks>0:
               tmp_args["ichunk"]=ichunk
               tmp_args["nchunks"]=nchunks
            tmp_all_args_list.append(tmp_args)
      all_args_list = tmp_all_args_list

   # Filter all_args_list
   all_args_list_filtered = []
   for dict_ent in all_args_list:
      fail_filters = False
      for filter_arg in extra_args:
         if filter_arg in dict_ent:
            if extra_args[filter_arg] != dict_ent[filter_arg]:
               fail_filters = True
         elif filter_arg == "output_file": # Can happen if output_file==short_name by default without explicit specification.
            if extra_args[filter_arg] != dict_ent["short_name"]:
               fail_filters = True
         else:
            dict_ent[filter_arg] = extra_args[filter_arg]
      if fail_filters:
         continue
      all_args_list_filtered.append(dict_ent)


   run_cmds = []
   if len(all_args_list_filtered)>0:
      if not os.path.isdir(jobs_cfgdir):
         os.makedirs(jobs_cfgdir)

      tarfilename = "ftanalysis.tar"
      tarfilepath = jobs_cfgdir + '/' + tarfilename
      if args.pkgtar is None:
         os.system("createFTAnalysisTarball.sh {}; mv {} {}".format(extratarcmd, tarfilename, tarfilepath))
      else:
         os.system("cp -f {} {}".format(args.pkgtar, tarfilepath))

      for all_args in all_args_list_filtered:
         strargs = get_strargs_from_dict(all_args)
         ichunk=-1
         nchunks=-1
         for xx in all_args:
            if xx=="ichunk":
               ichunk=int(all_args[xx])
            elif xx=="nchunks":
               nchunks=int(all_args[xx])

         output_file = ""
         if "output_file" in all_args:
            output_file = all_args["output_file"]
         else:
            output_file = all_args["short_name"]

         jobdirapp=""
         if nchunks>0:
            jobdirapp="_{}_{}".format(ichunk,nchunks)
         jobdir="{}/job_{}_{}_{}{}".format(jobs_cfgdir,args.exe,all_args["period"],output_file,jobdirapp)

         if args.overwrite and os.path.isdir(jobdir):
            os.system("rm -rf {}".format(jobdir))
         if not os.path.isdir(jobdir+"/Logs"):
            os.makedirs(jobdir+"/Logs")
         os.system("ln -sf {} {}/".format(tarfilepath, jobdir))

         jobargs={
            "BATCHQUEUE" : args.queue,
            "BATCHSCRIPT" : batch_script,
            "BATCHSITES" : allowed_sites,
            "TARFILE" : tarfilename,
            "OUTDIR" : jobdir,
            "CONDORSITE" : condor_site,
            "CONDOROUTDIR" : condor_outdir,
            "OUTLOG" : "Logs/log_job",
            "ERRLOG" : "Logs/err_job",
            "REQMEM" : args.memory,
            "REQDISK" : args.disk_space,
            "REQNCPUS" : args.ncpus,
            "JOBFLAVOR" : args.flavor,
            "EXE" : args.exe,
            "EXEARGS" : strargs
         }
         runCmd=str(
            "--dry --batchqueue={BATCHQUEUE} --batchscript={BATCHSCRIPT} --batchsites={BATCHSITES} --tarfile={TARFILE}" \
            " --condorsite={CONDORSITE} --condoroutdir={CONDOROUTDIR}" \
            " --outdir={OUTDIR} --outlog={OUTLOG} --errlog={ERRLOG}" \
            " --required_memory={REQMEM} --required_disk={REQDISK} --required_ncpus={REQNCPUS} --job_flavor={JOBFLAVOR}" \
            " --exe={EXE} {EXEARGS}"
         ).format(**jobargs)

         run_cmds.append(runCmd)

   print("Executing {} jobs...".format(len(run_cmds)))
   [ run_single(runCmd) for runCmd in run_cmds ]


if __name__ == '__main__':
   parser = ArgumentParser()

   parser.add_argument("extras", help="csv files with sample specifications (mandatory) and extra arguments to the executable", nargs="+")
   parser.add_argument("--exe", type=str, help="Executable name", required=True)
   parser.add_argument("--tag", type=str, help="Job tag", required=True)
   parser.add_argument("--nthreads", type=int, default=1, help="Number of threads to use in order to configure the batch jobs. Default=1.", required=False)
   parser.add_argument("--nevts", type=int, default=-1, help="Number of events per job. Default is -1 for no event splitting.", required=False)
   parser.add_argument("--queue", type=str, default="vanilla", help="Condor job queue", required=False)
   parser.add_argument("--flavor", type=str, default="tomorrow", help="Condor job flavor (walltime). Default=tomorrow", required=False)
   parser.add_argument("--sites", type=str, default=None, help="Condor job sites. Default is determined by the host name.", required=False)
   parser.add_argument("--memory", type=str, default="2048M", help="Condor job memory", required=False)
   parser.add_argument("--disk_space", type=str, default="5G", help="Condor job disk space", required=False)
   parser.add_argument("--ncpus", type=int, default=1, help="Number of CPUs to request for the condor job", required=False)
   parser.add_argument("--extra_tar", type=str, default=None, help="An extra tar file to include for specific execution needs", required=False)
   parser.add_argument("--group_eras", action="store_true", help="Group data sets for real data into data eras", required=False)
   parser.add_argument("--pkgtar", type=str, default=None, help="Use the precompiled tarball", required=False)
   parser.add_argument("--overwrite", action='store_true', help="Overwrite job directories even if they are present", required=False)

   args = parser.parse_args()

   run(args)
