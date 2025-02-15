import subprocess

# Define the GROMACS command
command = "wsl gmx pdb2gmx -f data/1aki.pdb -o data/1aki_processed.gro -water spc"

# Run the command and pass the input ('6' in this case)
try:
    # Pipe the input into the process, simulating user input for the force field choice
    result = subprocess.run(command, input="6\n", capture_output=True, text=True, shell=True)
    
    # Check if the command was successful
    if result.returncode == 0:
        print("Command Output:")
        print(result.stdout)
    else:
        print("Command failed with return code:", result.returncode)
        print("Command Errors:")
        print(result.stderr)

except Exception as e:
    print(f"Error running command: {e}")
