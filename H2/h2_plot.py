

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants as K
import pandas as pd
from mpl_toolkits import mplot3d

# Latex style
plt.style.use('default')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=20)
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
# Set ticks on both sides
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True


fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
ax.plot([0, 0], [0, 0], [0, 1], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([0, 0], [1,1], [0, 1],linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([0, 0], [0, 1], [0, 0], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([0, 1], [0, 0], [0, 0], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([0, 0], [0, 1], [1, 1],linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([0, 1], [0, 0], [1, 1],linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([1, 1], [0, 0], [0, 1], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([0, 1], [1, 1], [0, 0], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([1, 1], [1, 1], [0, 1], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([1, 1], [0, 1], [1, 1], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([0, 1], [1, 1], [1, 1], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([1, 1], [1, 1], [0, 1], linestyle = 'dashed',color = 'chocolate', linewidth=0.5)
ax.plot([1,1],[0,1],[0,0],'k--',linestyle = 'dashed',color = 'chocolate', linewidth=0.5)

ax.plot([0,0.5],[0,0.5],[0,0.5],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0.5,1],[0.5,1],[0.5,1],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0,0.5],[1,0.5],[0,0.5],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0.5,1],[0.5,0],[0.5,1],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0,0.5],[0,0.5],[1,0.5],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0,0.5],[1,0.5],[1,0.5],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0.5,1],[0.5,1],[0.5,0],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0.5,1],[0.5,0],[0.5,0],linestyle = 'dashed',color = 'silver', linewidth=1)
ax.plot([0.5,1],[0.5,0],[0.5,1],linestyle = 'dashed',color = 'silver', linewidth=1)



ax.scatter(0, 0, 0, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(0, 0, 1, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(0, 1, 0, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(1, 0, 0, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(1, 1, 1, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(1, 1, 0, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(1, 0, 1, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(0, 1, 1, c='chocolate', marker='o', alpha=1, s=650)
ax.scatter(0, 1, 1, c='chocolate', marker='o', alpha=1, s=50,label='Cu in $a_\mathrm{sub}$')
ax.scatter(0.5,0.5,0.5, c='silver', marker='o', alpha=1, s=50,label='Zn in $b_\mathrm{sub}$')
# ax.scatter(0.5,0.5,0.5, c='silver', marker='o', alpha=1, s=650)
# ax.scatter(1.5,0.5,0.5, c='silver', marker='o', alpha=1, s=650)
# ax.scatter(0.5,1.5,0.5, c='silver', marker='o', alpha=1, s=650)
# ax.scatter(0.5,0.5,1.5, c='silver', marker='o', alpha=1, s=650)
# ax.scatter(0.5,1.5,1.5, c='silver', marker='o', alpha=1, s=650)
# ax.scatter(1.5,0.5,1.5, c='silver', marker='o', alpha=1, s=650)
# ax.scatter(1.5,1.5,0.5, c='silver', marker='o', alpha=1, s=650)
# ax.scatter(1.5,1.5,1.5, c='silver', marker='o', alpha=1, s=650)


ax.padding = 10
ax.legend(loc='upper right',fontsize=13)



# cu_points = np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 1, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]])
# for i in range(len(cu_points)):
#     for j in range(i + 1, len(cu_points)):
#         ax.plot([cu_points[i, 0], cu_points[j, 0]], [cu_points[i, 1], cu_points[j, 1]], [cu_points[i, 2], cu_points[j, 2]], 'k-', linewidth=0.5)


ax.set_axis_off()
#plt.savefig('BCC_NN.pdf')#, bbox_inches='tight')  
plt.show()
