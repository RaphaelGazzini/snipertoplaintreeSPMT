import os 
import argparse
##################
#Argument for the job launcher 
# Create ArgumentParser object
parser = argparse.ArgumentParser(description="SniperToPlainTree jobLauncher ")

# Add arguments
parser.add_argument("--jobname", type=str,default="GetAndConv_", help="name of the job") 
parser.add_argument("--list", type=str, help="list of file to download with path", required=True)
parser.add_argument("--command", type=str, help="scp for IHEP cluster or xrdcp for eos",default ="scp")
parser.add_argument("--outpath", type=str,default="/sps/juno/llabit/", help="scp for IHEP cluster or xrdcp for eos") 
#parser.add_argument("--convert", action="store_true", help="convert to root")
parser.add_argument("--simcon", type=int, help="max simultaneous connections", default=10)
args = parser.parse_args()
#SBATCH OPTION: 
MEMORY = "5G"
NBTASK =1
NBCPU = 1 
JOBTIME = "5:00:00" # Hours:Minutes:Seconds
JOBNAME = args.jobname 
outpath = args.outpath
if not os.path.exists(outpath):
    os.makedirs(outpath, exist_ok=True)  # Creates the directory (including parent directories if needed)

listname = args.list 
listcorrname = listname+"_corr"
with open(listname,"r") as file1:
	num_lines = sum(1 for _ in file1)
print(num_lines)
with open(listcorrname,"w") as file2:
	with open(listname,"r") as file1:
		for l in file1:
			print(l)
			l2=l 
			l2 = l2.replace("esd/J25.4.3.a", "rtraw/")
			l2 = l2.replace("_J25.4.3.a.esd", ".rtraw")
			print(l2) 
			file2.write(l2)
#count=1
#for line in file1:
#	print(line)
JobName = JOBNAME#+str(count)
logname = JobName+"%a.log"
os.system("touch "+JobName+".sh")
f = open(JobName+".sh", "w")
f.write("#!/bin/sh \n")
f.write("# SLURM options:\n")
#	f.write("#SBATCH --requeue\n")
f.write("#SBATCH --job-name="+JobName+"\n")
f.write("#SBATCH --time "+str(JOBTIME)+"\n")
f.write("#SBATCH --cpus-per-task "+str(NBCPU)+"\n")
f.write("#SBATCH --mem "+str(MEMORY)+"\n")
f.write("#SBATCH --partition htc\n")
f.write("#SBATCH --output="+logname+"\n")
f.write("#SBATCH --array=1-"+str(num_lines)+"%"+str(args.simcon)+"\n")
f.write("#SBATCH --mail-user=loic.labit@iphc.cnrs.fr\n")
f.write("#SBATCH --mail-type=FAIL\n")   # options: BEGIN, END, FAIL, REQUEUE, ALL
f.write("# Commands to be submitted:\n")
f.write("export X509_USER_PROXY=$HOME/myproxy\n") 
f.write("line=$(sed -n ${SLURM_ARRAY_TASK_ID}p "+listname+")\n") 
f.write("linecorr=$(sed -n ${SLURM_ARRAY_TASK_ID}p "+listcorrname+")\n") 
f.write("filename=$(basename $linecorr | sed 's|.rtraw|.root|')\n")
f.write("name=${filename%.*}\n") 
command ="python run.py --input root://junoeos01.ihep.ac.cn/"+"$line"+" --input-correlations root://junoeos01.ihep.ac.cn/"+"$linecorr"+" --output "+outpath+"/$filename"
f.write(command)
f.close()
os.system("sleep 0.5")
os.system("sbatch "+JobName+".sh")
#	count=count+1
