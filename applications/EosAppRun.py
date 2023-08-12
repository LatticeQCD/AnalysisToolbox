import subprocess

def run_command(command):
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError:
        print(f"Failed to execute command: {command}")

def main():
    command1 = 'python main_HotQCDEoS.py --EosType "fixedsnB" --snB 50'
    command2 = 'python main_HotQCDEoS.py --EosType "fixedmuB" --muBdivT 2.0'
    
    run_command(command1)
    run_command(command2)

if __name__ == "__main__":
    main()

