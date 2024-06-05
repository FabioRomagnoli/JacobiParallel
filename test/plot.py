import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Plot performance data.')
parser.add_argument('filename', type=str, help='The CSV file containing the performance data')
args = parser.parse_args()

# Load data
try:
    data = pd.read_csv("data/" + args.filename + '.csv')
except Exception as e:
    print(f"Error reading the file: {e}")
    exit(1)

# Check if the required columns exist
required_columns = {'n_cores', 'threads', 'n_grid_points', 'time_elapsed', 'error'}
missing_columns = required_columns - set(data.columns)
if missing_columns:
    print(f"Missing columns in the data file: {missing_columns}")
    exit(1)

# Define colors for different numbers of cores
colors = {
    1: 'red',
    2: 'green',
    4: 'blue'
}

linestyles = {
    1: 'solid',
    2: 'dashed',
    4: 'dotted'
}


# Plotting
fig, ax = plt.subplots()

for n_cores in data['n_cores'].unique():
    for thread in data['threads'].unique():
        subset = data[(data['n_cores'] == n_cores) & (data['threads'] == thread)]
        ax.plot(subset['n_grid_points'], subset['time_elapsed'], label=f'{n_cores} cores, {thread} threads', 
                linestyle=linestyles.get(thread, 'solid'), color=colors.get(n_cores, 'black'))

ax.set_xscale('log', base=2)
ax.set_xlabel('Number of Grid Points')
ax.set_ylabel('Time Elapsed (s)')
ax.set_title('Performance Analysis')
ax.legend()
ax.grid(True)

plt.savefig('plots/performance_' + args.filename + '.png')
print(f"Plot saved: plots/performance_" + args.filename + ".png")

