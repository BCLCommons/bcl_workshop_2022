###############
#
# Simple plot script
#
#############

import sys
import numpy as np
import matplotlib.pyplot as plt

var1 = np.loadtxt(sys.argv[1])
var2 = np.loadtxt(sys.argv[2])
var3 = np.loadtxt(sys.argv[3])
output = sys.argv[4]
x = np.arange(len(var1))

plt.scatter(x, var1, c="blue", edgecolors="black",s=15, marker="o",alpha=0.50)
plt.scatter(x, var2, c="gray", edgecolors="black",s=15, marker="s",alpha=0.50)
plt.scatter(x, var3, c="green", edgecolors="black",s=15, marker="*",alpha=0.50)
plt.xlabel("Bin #")
plt.ylabel("Variance")
plt.title("3daSmoothSign",fontsize=10)
plt.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)
plt.tight_layout
plt.legend(("Dasatinib","Amprenavir","Ethinyl Estradiol"))
plt.savefig(output,dpi=300)
