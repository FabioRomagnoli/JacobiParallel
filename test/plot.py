import pandas as pd
import matplotlib.pyplot as plt

# Load data
try:
    data = pd.read_csv('./test/data/output.csv')
    print("Columns in DataFrame:", data.columns)
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

# Plotting
fig, ax = plt.subplots()

for n_cores in data['n_cores'].unique():
    for thread in data['threads'].unique():
        subset = data[(data['n_cores'] == n_cores) & (data['threads'] == thread)]
        linestyle = '-' if thread == 1 else '--'
        ax.plot(subset['n_grid_points'], subset['time_elapsed'], label=f'{n_cores} cores, {thread} threads', 
                linestyle=linestyle, color=colors.get(n_cores, 'black'))

ax.set_xscale('log', base=2)
ax.set_xlabel('Number of Grid Points')
ax.set_ylabel('Time Elapsed (s)')
ax.set_title('Performance Analysis')
ax.legend()
ax.grid(True)

plt.savefig('performance_plot.png')
plt.show()
