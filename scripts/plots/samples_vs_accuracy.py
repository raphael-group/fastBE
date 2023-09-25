import json
import os
import matplotlib.pyplot as plt
import seaborn as sns

# Initialize empty lists to store the numbers and objective values
numbers = []
objective_values = []

# Path to the directory where JSON files are stored
directory_path = "examples/"

# Loop through all files in the directory
for filename in os.listdir(directory_path):
    # Filter out files that don't match the pattern "alleleminima_*_results.json"
    if not filename.startswith("alleleminima_G_weighted") or not filename.endswith("_results.json"):
        continue

    # Extract the number from the filename
    number = int(filename.split("_")[3])

    # Open and read the JSON file
    with open(os.path.join(directory_path, filename), "r") as f:
        data = json.load(f)

    # Extract the objective value from the JSON data
    objective_value = data.get("objective_value", None)

    # Append the number and objective value to the lists
    if objective_value is not None:
        numbers.append(number)
        objective_values.append(objective_value)

# Create a Seaborn plot
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.scatterplot(x=numbers, y=objective_values)
sns.lineplot(x=numbers, y=objective_values)

# Add labels and title
plt.xlabel("Number of Random Samples")
plt.ylabel("$L_1(F, B)$")

# Show the plot
plt.show()

