import csv
import math
import numpy as np

# File path
csv_path = "convergence_results.csv"

# Read CSV data
N_values, error_p1, error_p2 = [], [], []

with open(csv_path, newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip header
    for row in reader:
        N_values.append(int(row[0]))
        error_p1.append(float(row[1]))
        error_p2.append(float(row[2]))

# Compute characteristic mesh size h = 1 / N
h_values = [1 / N for N in N_values]

# Compute convergence error rate

rates_p1 = [0] + [math.log(error_p1[i] / error_p1[i + 1]) / (math.log(N_values[i + 1] / N_values[i]) / 2) for i in range(len(error_p1) - 1)]
rates_p2 = [0] + [math.log(error_p2[i] / error_p2[i + 1]) / (math.log(N_values[i + 1] / N_values[i]) / 2) for i in range(len(error_p2) - 1)]

rates_p1[0] = 1
rates_p2[0] = 1

# Print the table
print(f"{'N':<10}{'h':<15}{'error_p1':<20}{'error_p2':<20}{'rate_p1':<15}{'rate_p2':<15}")
print("=" * 125)
for i in range(len(N_values)):
    print(f"{N_values[i]:<10}{h_values[i]:<15.8f}{error_p1[i]:<20.8e}{error_p2[i]:<20.8e}{rates_p1[i]:<15.6f}{rates_p2[i]:<15.6f}")
