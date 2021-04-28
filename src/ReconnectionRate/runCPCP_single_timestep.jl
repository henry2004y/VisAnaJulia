# Control script for calculating CPCP.
#
# Hongyang Zhou, hyzhou@umich.edu 11/18/2019

include("fieldline.jl")
################################################################################
# Main

# Step 1: create seeding points
createSeeds("seeds.txt", "../data/")

# Step 2: field line tracing (in Paraview)
pvpython = `pvpython ../scripts/streamtrace_save.py`;
run(pvpython)

# Step 3: find magnetopause boundary
file_stream = "../data/streamline.txt";
x, y, z = findBoundary(file_stream, true);
plot(x,y);
saveBoundary(x, y, z, "../data/boundary.txt");

# Step 4: interpolate onto the boundary points (in Paraview)
pvpython = `pvpython ../scripts/Interpolate_On_Boundary_batch.py`;
run(pvpython)

# Step 5: 2D line integral
file_integration = "../data/boundary_value.txt";
ϕ = integrate_along_boundary(file_integration, true)

file_integration = "../data/boundary_value_hall.txt";
ϕ_hall= integrate_along_boundary_hall(file_integration, true)
