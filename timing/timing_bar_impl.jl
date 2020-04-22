# Timings from implicit 256 core runs

using PyPlot

# 1 thread, 4 thread, 16 thread
# advance, krylov, Jacobian, facevalue, faceflux, exchange_message
t1 = [147.05, 123.81, 17.85, 3.28, 24.91, 4.96] ./ 3
t2 = [153.08, 129.88, 17.50, 3.19, 23.99, 4.97] ./ 3
t3 = [198.51, 172.18, 17.51, 3.23, 24.04, 12.99] ./ 3

t1_percent = [100.0, 84.2, 12.1, 2.2, 17.0, 3.4]
t2_percent = [100.0, 84.9, 11.4, 2.1, 15.7, 3.2]
t3_percent = [100.0, 86.7, 8.9, 1.6, 12.1, 6.5]


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
a2 = bar(r2, t2, width=barWidth, edgecolor="white", label="4 threads * 64 MPI")
a3 = bar(r3, t3, width=barWidth, edgecolor="white", label="16 threads * 16 MPI")

# Add xticks on the middle of the group bars
xlabel("module", fontweight=:bold, fontsize=18)
xticks([r + barWidth for r in 0:length(t1)-1],
   ["Advance", "Krylov", "Jacobian", "FaceValue", "FaceFlux", "MessagePass"],
   fontsize=18)
plt.minorticks_on()
ax.tick_params(axis=:x, which=:minor, bottom=false)
ax.yaxis.grid(true, linestyle="--")
#ax.yaxis.grid(true, which="minor")
ylabel("time [s]", fontsize=18)
ax.set_title("Timings per time step for 256 cores, implicit 2nd Order Linde + BiCGSTAB",fontsize=18)

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
