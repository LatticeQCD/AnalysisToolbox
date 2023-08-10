import subprocess
from multiprocessing import Pool
import argparse
from multiprocessing import freeze_support


def task(temp):
    subprocess.run(["python3", "main_HRG_LCP.py", "--r", str(r), "--hadron_file", filepath, "--models", "QM", "--T", str(temp)])

parser = argparse.ArgumentParser(description='Run different modes of the script.')
parser.add_argument('runMode', type=int, choices=[0, 1, 2], help='The run mode of the script.')
args = parser.parse_args()

r = 0.4
filepath = "../latqcdtools/physics/HRGtables/QM_hadron_list_ext_strange_2020.txt"
NTASKS = 10
        
def main():
        if args.runMode == 0:
            temps = range(120, 165)
            with Pool(NTASKS) as pool:
                pool.map(task, temps)
        
        elif args.runMode == 1:
            subprocess.run(["python3", "main_HRG_measure.py", "--hadron_file", filepath, "--models", "QM", "--LCP_file", "HRG_LCP_T100.0_r0.4", "--bqsc", "2000", "--obs", "chi"])
        
        elif args.runMode == 2:
            for obs in ["energy", "specificheat"]:
                for mb in [0.0, 1.0, 1.5, 2.0, 2.5]:
                    subprocess.run(["python3", "main_strangeness_neutral.py", "--hadron_file", filepath, "--fixedmuBNszerofile", f"HRG_fixedmuBT{mb}_r0.4QMHRG2020_BI", "--obs", obs, "--r", str(r), "--tag", "QMHRG2020_BI_Nszero"])
        else:
            print("Invalid runMode")


if __name__ == '__main__':
    freeze_support()
    main()
