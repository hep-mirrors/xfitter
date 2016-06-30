#!/usr/bin/env python

import plot_lhapdf
import lhapdf
from numpy import *
from matplotlib.pyplot import *

# nicer plots
rcParams['ytick.labelsize'] = 14
rcParams['xtick.labelsize'] = 14
rc('font',family='serif')
#import seaborn
#seaborn.set_style(style="white")

# load default PDF
aa = plot_lhapdf.loadPDF()

figure()
x = 0.3
qq = arange(1.,100.,0.1)
plot(qq,plot_lhapdf.fq(aa,1,x,qq)-plot_lhapdf.fq(aa,-1,x,qq))
xlabel('q, GeV',size=20)
ylabel('$u_v$',size=20)

savefig("uv.pdf")

figure()
xx = arange(-4.,0.,0.02)
q = 1.9
plot(xx, plot_lhapdf.fx(aa,1,10**xx,q),label="$u$",linewidth=2)
plot(xx, plot_lhapdf.fx(aa,2,10**xx,q),label="$d$",linewidth=2)
plot(xx, plot_lhapdf.fx(aa,3,10**xx,q),label="$s$",linewidth=2)
axis([-4.,0.,0.,1.05])
legend(fontsize=18)
xlabel('$ \log_{10} x$', size=20)
savefig("pdfs.pdf")

