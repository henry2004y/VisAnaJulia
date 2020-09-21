using PyPlot, PyCall

np = pyimport("numpy")

origin = "lower"

delta = 0.025

x = y = np.arange(-3.0, 3.01, delta)
X, Y = np.meshgrid(x, y)
Z1 = @. np.exp(-X^2 - Y^2)
Z2 = @. np.exp(-(X - 1)^2 - (Y - 1)^2)
Z = (Z1 .- Z2) .* 2

nr, nc = size(Z)

# mask a circle in the middle:
interior = @. sqrt(X^2 + Y^2) < 0.5
Z[interior] .= NaN

# We are using automatic selection of contour levels;
# this is usually not such a good idea, because they don't
# occur on nice boundaries, but we do it here for purposes
# of illustration.

fig1, ax2 = plt.subplots(constrained_layout=true)
CS = ax2.contourf(X, Y, Z, 10, cmap=plt.cm.bone, origin=origin)

# Note that in the following, we explicitly pass in a subset of
# the contour levels used for the filled contours.  Alternatively,
# We could pass in additional levels to provide extra resolution,
# or leave out the levels kwarg to use all of the original levels.

CS2 = ax2.contour(CS, levels=CS.levels[1:2:end], colors="r", origin=origin)

ax2.set_title("1 masked region")
ax2.set_xlabel("word length anomaly")
ax2.set_ylabel("sentence length anomaly")

# Make a colorbar for the ContourSet returned by the contourf call.
cbar = fig1.colorbar(CS)
cbar.ax.set_ylabel("verbosity coefficient")
# Add the contour line levels to the colorbar
cbar.add_lines(CS2)

plt.show()

