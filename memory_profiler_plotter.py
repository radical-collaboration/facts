import os
import matplotlib.pyplot as plt

def consolidate_logs(directory): 
    paths = []

    for root, dirs, files in os.walk(directory):
        for file in files:
            if "mem" in file: 
                # print(os.path.join(root, file))
                paths.append(os.path.join(root, file))

    paths.sort()

    memory_usage = {}
    mem_usages = {"fair": [], "Grls": [], "AIS": [], "stereo": [], "ESL": []}

    for path in paths:
        module = "fair"
        for key in mem_usages.keys(): 
            if key in path: 
                module = key

        with open(path, 'r') as file: 
            file.readline()      # skip first line
            line = file.readline()

            while line: 
                # memory_usage.append(float(line.split()[1]))
                timestamp = round(float(line.split()[2]), 1)
                memory = float(line.split()[1])

                mem_usages[module].append(memory)

                memory_usage[timestamp] = memory + memory_usage.get(timestamp, 0)

                line = file.readline()

        for key in mem_usages.keys():
            mem_usages[key].sort()

    print(f"Max Memory Usages: {[(key, mem_usages[key][-1]) for key in mem_usages.keys()]}")
    print(f"95th Percentile: {[(key, mem_usages[key][int(0.95 * (len(mem_usages[key]) - 1))]) for key in mem_usages.keys()]}")
                
    plt.figure()
    plt.plot(list(range(len(memory_usage.keys()))), [memory for timestamp, memory in sorted(memory_usage.items())])
    plt.title("Memory Usage Over Time of 1k Locations and 2k Samples")
    plt.xlabel("Time (Tenths of Second)")
    plt.ylabel("Megabytes")
    plt.savefig("memory_usage.png")


if __name__ == "__main__": 
    # consolidate_logs("~/radical.pilot.sandbox/full_12Locs_2kSamples")
    consolidate_logs(os.path.join(os.path.expanduser("~"), "radical.pilot.sandbox/full_1kLocs_2kSamples"))
    # consolidate_logs(os.path.join(os.path.expanduser("~"), "radical.pilot.sandbox/fairOnly_1Locations_5kSamples"))