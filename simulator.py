import requests

ligand_name = "your_ligand_name"  # Replace with actual ligand name
response = requests.post("http://localhost:5000/update", json={"ligand_name": ligand_name})
print(response.json())  # Should print {"message": "Output updated"}
