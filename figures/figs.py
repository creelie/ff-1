import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
plt.rcParams.update({"font.family":"serif","font.size":11,"mathtext.fontset":"cm"})

# Figure 1: the subschemes S_{q0}, S_{q1}, S_H inside P^1 over Spec Z (schematic, 3D)
fig=plt.figure(figsize=(7.2,4.4)); ax=fig.add_subplot(111,projection="3d")
u=np.linspace(0,8.4,300)
for p in range(0,9):
    t=np.linspace(-np.pi/2,np.pi/2,80)
    ax.plot(p+0*t,np.cos(t),np.sin(t),color="0.80",lw=0.8)
def band(ang,col,lab=None,lw=2.2):
    ax.plot(u,np.cos(ang),np.sin(ang),color=col,lw=lw,label=lab)
b0=-1.05+0.06*np.sin(0.8*u); b0b=-0.55+0.06*np.cos(0.6*u)
band(b0,"#1f4e9c",r"$S_{q_0}$"); band(b0b,"#1f4e9c")
bump=1-np.exp(-((u-2.0)/0.9)**2)
band(b0+0.17*bump,"#b5402a",r"$S_{q_1}$"); band(b0b-0.15-0.02*np.sin(0.7*u),"#b5402a")
for j,a in enumerate([0.25,0.55,0.85,1.15]):
    band(a+0.05*np.sin(0.5*u+j),"#2f7d32",r"$S_H$" if j==0 else None,lw=1.7)
ax.scatter([2.0],[np.cos(b0[np.argmin(abs(u-2.0))])],[np.sin(b0[np.argmin(abs(u-2.0))])],color="k",s=22,zorder=5)
ax.text(2.0,np.cos(-1.05)+0.05,np.sin(-1.05)-0.32,r"$T$",fontsize=12)
ax.set_xlabel(r"$\mathrm{Spec}\,\mathbb{Z}$",labelpad=8); ax.set_yticks([]); ax.set_zticks([])
ax.set_xticks(range(0,9)); ax.set_xticklabels(["","2","3","5","7","11","13","17","19"])
ax.view_init(elev=8,azim=-75); ax.set_box_aspect((2.8,0.8,1.2))
ax.legend(loc="upper left",frameon=False)
plt.tight_layout(); plt.savefig("fig_scheme.png",dpi=220,bbox_inches="tight",pad_inches=0.05); plt.close()

# Figure 2: degrees covered
N=120
fig,ax=plt.subplots(figsize=(7.2,2.9))
cols=[]
for n in range(4,N+1):
    if n==4: c="#888888"
    elif n%2==0: c="#1f4e9c"
    elif (n-2)%3==0: c="#b5402a"
    elif (n-2)%5==0: c="#e08a1e"
    elif (n-2)%7==0: c="#2f7d32"
    else: c="white"
    i=n-4; x=i%20; y=-(i//20)
    ax.add_patch(plt.Rectangle((x,y),0.92,0.92,facecolor=c,edgecolor="0.4",lw=0.6))
    ax.text(x+0.46,y+0.46,str(n),ha="center",va="center",fontsize=7,color="white" if c!="white" else "black")
ax.set_xlim(-0.2,20.2); ax.set_ylim(-5.2,1.1); ax.axis("off")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor="#888888",label="BEGRS"),Patch(facecolor="#1f4e9c",label=r"$n$ even"),
  Patch(facecolor="#b5402a",label=r"$3\mid n-2$"),Patch(facecolor="#e08a1e",label=r"$5\mid n-2$"),
  Patch(facecolor="#2f7d32",label=r"$7\mid n-2$"),Patch(facecolor="white",edgecolor="0.4",label="open")],
  loc="upper center",bbox_to_anchor=(0.5,1.28),ncol=6,frameon=False,fontsize=8.5)
plt.tight_layout(); plt.savefig("fig_degrees.png",dpi=220,bbox_inches="tight",pad_inches=0.05); plt.close()
