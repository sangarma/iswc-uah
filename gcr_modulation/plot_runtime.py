import matplotlib.pyplot as plt
import numpy as np
cores = [1,2,4,8]
runtime = [25,14,8,5]
plt.plot(cores,runtime, label = "runtime")
plt.plot(cores,25/np.array(cores), label = "1/n")
plt.xlabel("number of cores")
plt.ylabel("runtime")
plt.legend()
plt.show()