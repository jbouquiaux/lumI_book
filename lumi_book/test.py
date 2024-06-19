from abipy.lumi.deltaSCF import DeltaSCF
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import numpy as np


files=["lumi/flow_deltaSCF/w0/t2/outdata/out_GSR.nc",
                                        "lumi/flow_deltaSCF/w0/t3/outdata/out_GSR.nc",
                                        "lumi/flow_deltaSCF/w0/t4/outdata/out_GSR.nc",
                                        "lumi/flow_deltaSCF/w0/t5/outdata/out_GSR.nc",]
results=DeltaSCF.from_four_points_file(files) 
stru=results.structure_ex()

fig,axs=plt.subplots(2,1)

stru.plot_xrd(ax=axs[0]);
stru.plot_xrd(ax=axs[1]);




