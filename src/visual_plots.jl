# Using user recipes from Plots.

using RecipesBase

#@userplot mycontour

# Build a recipe which acts on a custom type.
# The function name here is meaningless: it is only used to process a unique 
# set of types early in the pipeline.
# It will work on all functions given the correct dimensions, e.g.
# plot(data, "p")
# contourf(data, "Mx")
@recipe function f(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1)

   ndim = data.head.ndim

   if ndim == 1
      VarIndex_ = findindex(data, var)
      x, w = data.x, data.w
      y = w[:,VarIndex_]

      @series begin
         seriestype --> :path
         x, y
      end
   elseif ndim == 2
      xi, yi, wi = getdata(data, var, plotrange, plotinterval)

      @series begin
         seriestype --> :contourf  # use := if you want to force it
         xi, yi, wi
      end
   end
end
