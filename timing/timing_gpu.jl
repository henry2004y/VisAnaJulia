# Timings for the skeleton BATSRUS for porting to GPU with OpenACC.
# 4 subplots, each of block size 8^3, 16^3, 32^3, and 64^3.
# Each subplot contains a bar plot, with nBlock as x axis, timing as y axis,
# and different columns represent nWork.
# Column have two colors, blue represent CPU, and orange represent GPU.
#
# Hongyang Zhou, hyzhou@umich.edu 04/22/2020

using PyPlot, CSV, DataFrames

#
df = DataFrame(CSV.File("results/timing_gpu.txt"))
t1 = df[df.nCell.==8,:]
t2 = df[df.nCell.==16,:]
t3 = df[df.nCell.==32,:]
t4 = df[df.nCell.==64,:]

t1_cpu = t1[t1.type.=="CPU",:]

# set width of bar
const barWidth = 0.3

# Set position of bar on X axis
r1 = collect(0:size(t1,1)-1)
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

fig, ax = subplots(2,2,figsize=(14,6))
#axis("tight")

# Make the plot
#a1 = ax[1,1].bar(r1, t1[:,:time], width=barWidth, edgecolor="white", label="1 thread * 256 MPI")
#a2 = ax[1,2].bar(r2, t2[:,:time], width=barWidth, edgecolor="white", label="8 threads * 32 MPI")
#a3 = ax[2,1].bar(r3, t3[:,:time], width=barWidth, edgecolor="white", label="32 threads * 8 MPI")
#a3 = ax[2,2].bar(r3, t4[:,:time], width=barWidth, edgecolor="white", label="32 threads * 8 MPI")
#=
# Add xticks on the middle of the group bars
xlabel("module", fontweight=:bold, fontsize=18)
xticks([r + barWidth for r in 0:length(t1)-1],
   ["Advance", "FaceValue", "FaceFlux", "Update", "UpdateCheck", "MessagePass"],
   fontsize=18)
plt.minorticks_on()
ax.tick_params(axis=:x, which=:minor, bottom=false)
ax.yaxis.grid(true, linestyle="--")
#ax.yaxis.grid(true, which="minor")
ylabel("time [s]", fontsize=18)
ax.set_title("Timings per time step for 256 cores, explicit 2nd Order Linde + mc3 limiter",fontsize=18)

# Create legend & Show graphic
legend(fontsize=18)

function autolabel(rects, t_percent)
   """Attach a text label above each bar in *rects*, displaying its height."""
   for (i, rect) in enumerate(rects)
      i == 1 && continue
      height = rect.get_height()
      percent = t_percent[i]
      ax.annotate(string(percent),
         xy=(rect.get_x() + rect.get_width() / 2, height),
         xytext=(0, 3),  # 3 points vertical offset
         textcoords="offset points",
         ha=:center, va=:bottom)
   end
end

autolabel(a1, t1_percent)
autolabel(a2, t2_percent)
autolabel(a3, t3_percent)
=#
