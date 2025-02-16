import os

def load_pdb_file(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    with open(file_path, 'r') as file:
        return file.read()

def save_results(file_path, data):
    with open(file_path, 'w') as file:
        file.write(data)

def list_files_in_directory(directory):
    if not os.path.exists(directory):
        raise FileNotFoundError(f"The directory {directory} does not exist.")
    return [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)