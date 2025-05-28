using Literate
using Base
execute = true
# ENV["JULIA_DEBUG"]="Literate"
Literate.markdown("./extrapolation.jl", pwd()*"/../docs/src/examples/extrapolation/"; execute)
Literate.notebook("./extrapolation.jl", pwd()*"/../examples/"; execute)
Literate.markdown("./birkhoff_averaging.jl", pwd()*"/../docs/src/examples/birkhoff_averaging/"; execute)
Literate.notebook("./birkhoff_averaging.jl", pwd()*"/../examples/"; execute)
Literate.markdown("./kernel.jl", pwd()*"/../docs/src/examples/kernel/"; execute)
Literate.notebook("./kernel.jl", pwd()*"/../examples/"; execute)
if execute
	Base.rm("./3D_ER3BP/Trojan.jld2")
end
Literate.markdown("./3d_extrapolation.jl", pwd()*"/../docs/src/examples/3d_extrapolation/"; execute)
Literate.notebook("./3d_extrapolation.jl", pwd()*"/../examples/"; execute)
