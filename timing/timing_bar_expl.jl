# Timings from Strong scaling 256 core runs

using PyPlot

# 1 thread, 8 thread, 32 thread
t1 = [6.62, 1.3, 4.0, 0.24, 0.21, 0.69] ./ 5
t2 = [6.85, 1.29, 4.01, 0.24, 0.2, 0.81] ./ 5
t3 = [10.69, 1.34, 4.05, 0.27, 0.88, 3.87] ./ 5

t1_percent = [99.34, 19.6, 60.4, 3.7, 3.1, 10.5]
t2_percent = [99.11, 18.9, 58.5, 3.5, 3.0, 11.8]
t3_percent = [99.53, 12.5, 37.9, 2.5, 8.2, 36.18]

#=
t_all = [6.62, 6.85, 10.69]
t_facevalue = [1.3, 1.29, 1.34]
t_faceflux = [4.0, 4.01, 4.05]
t_source = [0.07, 0.07, 0.08]
t_update = [0.24, 0.24, 0.27]
t_updatecheck = [0.21, 0.2, 0.88]
t_messagepass = [0.69, 0.81, 3.87]
=#

# set width of bar
const barWidth = 0.3

# Set position of bar on X axis
r1 = collect(0:length(t1)-1)
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]

fig = figure("bar_timing",figsize=(14,6))
ax = plt.axes()
#axis("tight")

# Make the plot
a1 = bar(r1, t1, width=barWidth, edgecolor="white", label="1 thread * 256 MPI")
a2 = bar(r2, t2, width=barWidth, edgecolor="white", label="8 threads * 32 MPI")
a3 = bar(r3, t3, width=barWidth, edgecolor="white", label="32 threads * 8 MPI")

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
