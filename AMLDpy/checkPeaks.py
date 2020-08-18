mFile = pd.read_csv('/Users/emilywilliams/Documents/GitHub/AMLD_Driving_Data/coDrive_new/coDrive_10pc_45_102_median/FinalShpFiles/mainThing.csv')

import matplotlib.pyplot as plt

# An "interface" to matplotlib.axes.Axes.hist() method
n, bins, patches = plt.hist(x=mFile.OP_DISTANCE, bins=20, color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('OP Distance (m)')
plt.ylabel('Frequency')
plt.title('Observed Peak Distance - co drives')
plt.text(23, 45, r'$Threshold=10pc, Baseline=median$')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
plt.show()

longPk = mFile.loc[mFile.OP_DISTANCE > 200,:]


fig1, ax1 = plt.subplots()
ax1.set_title('Basic Plot')
ax1.boxplot(mFile.OP_DISTANCE)

fig2, ax2 = plt.subplots()
ax2.set_title('Notched boxes')
ax2.boxplot(float(mFile.OP_DISTANCE), notch=True)

from functools import reduce
test = [1,20,2,23,93]
reduce(lambda x,y: x+y,test)