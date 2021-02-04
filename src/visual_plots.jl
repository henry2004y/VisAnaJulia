# Using user recipes from Plots.

using RecipesBase
using UnitfulBatsrus, UnitfulRecipes
# UnitfulBatsrus is an unregistered package!

#@userplot mycontour

# Build a recipe which acts on a custom type.
# The function name here is meaningless: it is only used to process a unique 
# set of types early in the pipeline.
# It will work on all functions given the correct dimensions, e.g.
# plot(data, "p")
# contourf(data, "Mx", xlabel="x")
@recipe function f(data::Data, var::AbstractString;
   plotrange=[-Inf,Inf,-Inf,Inf], plotinterval=0.1)

   ndim = data.head.ndim

   if startswith(data.head.headline, "normalized")
      hasunits = false
   else
      hasunits = true
      unitw = getunit(data, var)
   end

   if ndim == 1
      VarIndex_ = findindex(data, var)
      if hasunits
         unitx = getunit(data, data.head.variables[1])
         x = data.x .* unitx
         w = data.w
         y = w[:,VarIndex_] .* unitw
      else
         x, w = data.x, data.w
         y = w[:,VarIndex_]
      end

      @series begin
         seriestype --> :path
         x, y
      end
   elseif ndim == 2
      x, y, w = getdata(data, var, plotrange, plotinterval)
      unitx = getunit(data, data.head.variables[1])
      unity = getunit(data, data.head.variables[2])

      x *= unitx
      y *= unity
      w *= unitw

      @series begin
         seriestype --> :contourf  # use := if you want to force it
         x, y, w'
      end
   end
end


function getunit(data, var)

   # Batrus has a bug in the 2D cuts of 3D runs: it always outputs the 3
   # coordinate units in the headline. To temporarily deal with it, I plus
   # the index by 1.
   var_ = findfirst(x->x==lowercase(var), lowercase.(data.head.variables)) + 1
   isnothing(var_) && error("$(var) not found in file header variables!")
   var_unit_strs = split(data.head.headline)

   if var_unit_strs[var_] == "R"
      var_unit = bu"R"
   elseif var_unit_strs[var_] == "Mp/cc"
      var_unit = bu"amucc"
   elseif var_unit_strs[var_] == "uA/m2"
      var_unit = bu"ampm2"
   elseif var_unit_strs[var_] == "V/m2"
      var_unit = bu"vm2"
   else
      var_unit = UnitfulBatsrus.Unitful.uparse(var_unit_strs[var_])
   end
   var_unit
end

function getunits(data)

   var_unit_strs = split(data.head.headline)
   var_units = [] # needs to be improved!
   for var_unit_str in var_unit_strs
      if var_unit_str == "R"
         var_unit = bu"R"
      elseif var_unit_str == "Mp/cc"
         var_unit = bu"amucc"
      elseif var_unit_str == "uA/m2"
         var_unit = bu"ampm2"
      elseif var_unit_str == "V/m2"
         var_unit = bu"vm2"
      else
         var_unit = UnitfulBatsrus.Unitful.uparse(var_unit_str)
      end
      push!(var_units, var_unit)
   end
   var_units
end