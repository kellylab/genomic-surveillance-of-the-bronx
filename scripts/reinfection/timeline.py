import pandas as pd
import matplotlib.pyplot as plt

cts = pd.read_csv("data/external/ct_timeline.csv",index_col=0, header=None)
cts.columns = ['CT Values']
cts.index = pd.to_datetime(cts.index)   
cts.index.name = 'Date'
cts['Sampling'] = [True, False, False, False, False, False, True, False, False, True, False, False, False]   
dx = pd.DataFrame(index=pd.date_range(cts.index.min(),cts.index.max()))
dx['CT Values'] = cts['CT Values']
dx = dx.fillna(0.)
plt.clf()
fig, ax = plt.subplots()
cts.sort_index().plot(linewidth=2, legend=False,ax=ax, style='--')
cts.dropna().reset_index().plot.scatter(x='Date',y='CT Values', ax=ax, color='blue', s=80,linewidth=3)
(cts[cts['Sampling']]).reset_index().plot.scatter(x='Date',y='CT Values',ax=ax,color='green',s=130, linewidth=3)
# dx['CT Values'].plot.bar(ax=ax)
ax.set_ylabel("CT Values")
ax.set_xticks(cts.index) 
ax.set_xticklabels(cts.index.map(lambda x: f"{x.month}-{x.day}"),rotation=90)
ax.set_ylim([0,35])
plt.savefig("data/processed/reinfection/timeline.pdf")
plt.show()