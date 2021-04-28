using Glob, DelimitedFiles, PyPlot#, Plots


files = glob("boundary[0-9]*.txt", "data")

xy = [zeros(Float32, 0,2) for _ in 1:length(files)]

for i in 1:length(files)
   bc = readdlm(files[i], ',', Float32, header=true)
   xy[i] = bc[1][:,1:2]
end

figure()

for i ∈ 1:length(files)
   plot(xy[i][:,1],xy[i][:,2], label="boundary")
   xlim(-2.1,4.2)
   ylim(-3.0,3.0)
   title("t = $i")
   savefig("figure/PIC_test"*lpad(i,4,"0")*".png")
   clf()
end


#=
anim = @animate for i ∈ 1:length(files)
   plot(xy[i][:,1],xy[i][:,2], xlims=(-2.1,4.2), ylims=(-3.0,3.0), label="boundary")
   title!("t = $i")
end

gif(anim, "test.mp4", fps = 15)
=#
