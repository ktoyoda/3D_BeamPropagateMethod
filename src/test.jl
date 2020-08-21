using Plots
using KTOptical

um = 10^-6
KTOptical.setParam(5um, 0.0, 1.06um);

x = range(-10um,10um,step = 0.1um)
y = range(-10um,10um,step = 0.1um)
I = HG_I.(1,1,x,y')

typeof(x)
