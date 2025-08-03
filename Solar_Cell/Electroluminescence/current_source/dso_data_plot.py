# %%
import csv
import numpy as np
import matplotlib.pyplot as plt

# %%
file_path = "SDS2104X Plus_CSV_C1_1.csv"

time = []
y = []
with open(file_path, 'r') as f:
    data = csv.reader(f)
    for raw in data:
        time.append(raw[0])
        y.append(raw[1])

# %%
t = []
I = []
for i in range(12,len(time)):
    t.append(float(time[i]))
    I.append(float(y[i])) 

# %%
plt.plot(t,I)


