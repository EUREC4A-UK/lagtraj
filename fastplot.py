from pylab import *
import xarray as xr

my_scale=[-1e-4,1e-4]
my_norm = mpl.colors.Normalize(vmin=my_scale[0], vmax=my_scale[1])
levels = linspace(my_scale[0],my_scale[1],21)

ds=xr.open_dataset('ds_along_traj.nc')
subplot(2,1,1)
ds['dp_fdx'].transpose().plot.contourf(levels=levels,norm=my_norm,cmap='RdBu_r')
title('regression')
subplot(2,1,2)
ds['dp_fdx_bound'].transpose().plot.contourf(levels=levels,norm=my_norm,cmap='RdBu_r')
title('boundaries')
gcf().tight_layout()
show()
close()

my_scale_2=[-4.e-4,4.e-4]
my_norm_2 = mpl.colors.Normalize(vmin=my_scale_2[0],vmax=my_scale_2[1])
levels_2 = linspace(my_scale_2[0],my_scale_2[1],21)

subplot(2,1,1)
ds['dp_fdy'].transpose().plot.contourf(levels=levels_2,norm=my_norm_2,cmap='RdBu_r')
title('regression')
subplot(2,1,2)
ds['dp_fdy_bound'].transpose().plot.contourf(levels=levels_2,norm=my_norm_2,cmap='RdBu_r')
title('boundaries')
gcf().tight_layout()
show()
close()

my_scale_3=[-10.,10.]
my_norm_3 = mpl.colors.Normalize(vmin=my_scale_3[0],vmax=my_scale_3[1])
levels_3 = linspace(my_scale_3[0],my_scale_3[1],21)

subplot(2,1,1)
ds['u'].transpose().plot.contourf(levels=levels_3,norm=my_norm_3,cmap='RdBu_r')
title('point')
subplot(2,1,2)
ds['u_mean'].transpose().plot.contourf(levels=levels_3,norm=my_norm_3,cmap='RdBu_r')
title('mean')
gcf().tight_layout()
show()
