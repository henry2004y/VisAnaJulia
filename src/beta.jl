# Plot plasma β from PIC outputs.
#
# Hongyang Zhou, hyzhou@umich.edu

using Batsrus, VisAna, PyPlot, Printf

filename = "3d_box.out"
filehead, data, filelist = readdata(filename)

cutPlaneIndex = 65
VarIndex_ = 18

X = @view data.x[:,:,:,1]
Y = @view data.x[:,:,:,2]
Z = @view data.x[:,:,:,3]

fig, ax = plt.subplots(1, 1)
fig.set_size_inches(3.3,6)

W  = @view data.w[:,:,:,VarIndex_]
Bx = @view data.w[:,:,:,5]
By = @view data.w[:,:,:,6]
Bz = @view data.w[:,:,:,7]

cut1 = @view X[:,cutPlaneIndex,:]
cut2 = @view Z[:,cutPlaneIndex,:]
W    = @view W[:,cutPlaneIndex,:]
Bx   = @view Bx[:,cutPlaneIndex,:]
By   = @view By[:,cutPlaneIndex,:]
Bz   = @view Bz[:,cutPlaneIndex,:]
PB   = hypot.(Bx, By, Bz)

c = ax.contourf(cut1,cut2,W./PB)
fig.colorbar(c,ax=ax)
ax.axis("scaled")
title(filehead[:wnames][VarIndex_])

xlabel("x"); ylabel("z")

dim = [0.125, 0.013, 0.2, 0.045]
str = @sprintf "it=%d, time=%4.2f" filehead[:it] filehead[:time]
at = matplotlib.offsetbox.AnchoredText(str,
   loc="lower left", prop=Dict("size"=>8), frameon=true,
   bbox_to_anchor=(0., 1.),
   bbox_transform=ax.transAxes)
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(at)
