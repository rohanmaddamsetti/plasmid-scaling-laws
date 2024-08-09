import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Load the combined CSV file
combined_file_path = '../results/Combined.csv'
df_combined = pd.read_csv(combined_file_path)

# Group the replicon_length and CopyNumber by strain identity and compute total CopyNumber for each strain
grouped = df_combined.groupby('Strain')
total_copy_number = grouped['CopyNumber'].sum()

# Compute the average length weighted by the copyNumber for each strain
weighted_average_length = grouped.apply(lambda x: np.average(x['replicon_length'], weights=x['CopyNumber']))

# Generate log-log plot of total CopyNumber against weighted average length
fig, ax = plt.subplots()
ax.loglog(weighted_average_length, total_copy_number, 'o', label='Data')

# Perform linear regression on log-transformed data
slope, intercept, r_value, p_value, std_err = linregress(np.log(weighted_average_length), np.log(total_copy_number))
equation = f'y = {np.exp(intercept):.2f} * x^{slope:.2f}'
r_squared = r_value ** 2

# Plot the linear fit on log-log scale
x_fit = np.linspace(min(np.log(weighted_average_length)), max(np.log(weighted_average_length)), 100)
y_fit = np.exp(intercept) * np.exp(slope * x_fit)
ax.loglog(np.exp(x_fit), y_fit, 'r-', label='Linear Fit')

# Set plot labels and legend
ax.set_xlabel('Weighted Average Length (log scale)')
ax.set_ylabel('Total Copy Number (log scale)')
ax.set_title('Total Copy Number vs. Weighted Average Length')
ax.legend()

# Display the equation and R-squared value
text = f'Equation: {equation}\nR-squared: {r_squared:.2f}'
ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=10, verticalalignment='top')

# Show the plot
plt.show()
