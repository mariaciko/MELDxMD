import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


trace=np.loadtxt("trace.dat")
np.shape(trace)

trace = trace.transpose()
np.shape(trace)

print(trace[15])


# TEST WITH FAKE DATA

fake_data = [29, 29, 29, 0, 0, 0, 5, 7, 8, 0, 1, 28, 29, 0, 3, 2, 9, 0, 27, 0, 29, 0]

def check_if_29_or_0(runner):
    meta_traj = []
    for i in runner:
        if i== 0:
            meta_traj.append(i)
        if i == 29:
            meta_traj.append(i)
    return meta_traj

def clean_traj(meta_traj):
    clean_traj = []
    clean_traj.append(meta_traj[0])
    for i in range(1, len(meta_traj)):
        if meta_traj[i] != meta_traj[i-1]:
            clean_traj.append(meta_traj[i])            
    return clean_traj


first_output = check_if_29_or_0(fake_data)
first_output

cleaned_traj = clean_traj(first_output)
cleaned_traj

number_of_roundtrips = int((len(cleaned_traj) - 3)/2) +1
number_of_roundtrips



# ACTUAL DATA

# Pass actual data
trace=np.loadtxt("trace.dat")
trace.shape

trace = trace.transpose()
trace.shape
print(trace)

first_output = check_if_29_or_0(trace[0])
print("First", first_output)
cleaned_traj = clean_traj(first_output)
print("Second", cleaned_traj)
number_of_roundtrips = int((len(cleaned_traj) - 3)/2) +1
print("Number of roundtrips", number_of_roundtrips)

total_roundtrips = []
for i in trace:
    first_output = check_if_29_or_0(i)
    cleaned_traj = clean_traj(first_output)
    number_of_roundtrips = int((len(cleaned_traj) - 3)/2) +1
    print("Number of roundtrips", number_of_roundtrips)
    total_roundtrips.append(number_of_roundtrips)
    print(f'Runner {i}: Number of roundtrips = {number_of_roundtrips}')
print("Total round trips:", total_roundtrips)


# Compute the mean and standard deviation to identify outliers
mean_roundtrips = np.mean(total_roundtrips)
std_roundtrips = np.std(total_roundtrips)

# Set a threshold for defining outliers (e.g., below one standard deviation from the mean)
outlier_threshold = mean_roundtrips - std_roundtrips

# Initialize the figure for the histogram
plt.figure(figsize=(10, 6))

# Generate a list of colors: highlight runners with fewer roundtrips than the threshold
colors = ['red' if rt < outlier_threshold else 'blue' for rt in total_roundtrips]

# Create a bar plot (histogram) with runners on the x-axis and round trips on the y-axis
bars = plt.bar(range(len(total_roundtrips)), total_roundtrips, color=colors)

# Add labels and title
plt.xlabel('Runner Index', fontsize=12)
plt.ylabel('Number of Roundtrips', fontsize=12)
plt.title('Roundtrips per Runner', fontsize=14)

# Annotate the bars with the number of roundtrips
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2.0, height, '%d' % int(height), ha='center', va='bottom')

# Show the plot
plt.tight_layout()
plt.savefig('roundtrips.png')
plt.close()