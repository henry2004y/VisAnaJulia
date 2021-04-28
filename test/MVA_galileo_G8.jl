# Example minimum variance analysis (MVA) from Galileo satellite Ganymede flyby
# magnetometer data.
# G8 inbound magnetopausee LMN for interval 15:51:40-15:52:20
# I cannot repeat the results of Figure 6 in [Kivelson 1998]
# -0.72 -0.55 0.43 from [Kivelson 1998]
# -0.942 -0.187 0.278 from elliptic empirical model
#
# Hongyang Zhou, hyzhou@umich.edu 01/10/2020

using CSV, DataFrames, CodecZlib, Dates, PyPlot, Mmap

"""
    MVA_analysis(filename::AbstractString, index_::UnitRange, DoPlot=false)

This function is currently only written for the specific input file.
"""
function MVA_analysis(filename::AbstractString, index_::UnitRange, DoPlot=false)
   # Load magnetometer data Bx, By, Bz
   if filename[end-1:end] == "gz"
      df = CSV.File(transcode(GzipDecompressor, Mmap.mmap(filename));
         header=2, delim=" ", ignorerepeated=true) |> DataFrame
   else
      df = CSV.File(filename; header=2, delim=" ", ignorerepeated=true)
   end

   # inbound crossing interval
   F = MVA(df.Bx[index_], df.By[index_], df.Bz[index_])

   t = DateTime.(df.yr, df.month, df.day, df.hr, df.min, floor.(Int, df.sec),
      floor.(Int, 1e3 .* (df.sec - floor.(df.sec))) )

   if DoPlot
      # Original Cartesian coordinates
      figure(figsize=(12,5))
      subplot(411)
      plot(t, df.Bx, label=L"B_x")
      legend()
      subplot(412)
      plot(t, df.By, label=L"B_y")
      legend()
      subplot(413)
      plot(t, df.Bz, label=L"B_z")
      legend()
      subplot(414)
      plot(hypot.(df.Bx, df.By, df.Bz), label=L"B")
      legend()
      tight_layout()

      # Coordinate transformation
      id_ = 2000:2900
      BL = sum(F.vectors[:,1]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
      BM = sum(F.vectors[:,2]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
      BN = sum(F.vectors[:,3]' .* [df.Bx[id_] df.By[id_] df.Bz[id_]];dims=2)
      figure(figsize=(12,5))
      subplot(411)
      plot(t[id_], BL, label=L"B_L")
      legend()
      subplot(412)
      plot(t[id_], BM, label=L"B_M")
      legend()
      subplot(413)
      plot(t[id_], BN, label=L"B_N")
      legend()
      subplot(414)
      plot(t[id_], hypot.(BL, BM, BN), label=L"B")
      legend()
      tight_layout()
   end

   F
end