pyplot() # Good for publication quality plots (requires PyPlot.jl)

seaborn = palette(:seaborn_deep)

# General plot settings
default(fontfamily="serif", legendfontsize=10)
default(titlefontsize=12, grid=false, framestyle=:box)

default(xticks=:auto, yticks=:auto, guidefont=("serif", 12))
default(xtickfont=("serif", 10), ytickfont=("serif", 10))
default(xlabelfont=("serif", 12), ylabelfont=("serif", 12))