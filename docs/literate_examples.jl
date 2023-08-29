using Literate
using Base
execute = true
# ENV["JULIA_DEBUG"]="Literate"
Literate.notebook("./extrapolation.jl", pwd()*"/../examples/"; execute)
Literate.markdown("./extrapolation.jl", pwd()*"/../docs/src/examples/extrapolation/"; execute)
Literate.notebook("./birkhoff_averaging.jl", pwd()*"/../examples/"; execute)
Literate.markdown("./birkhoff_averaging.jl", pwd()*"/../docs/src/examples/birkhoff_averaging/"; execute)
Literate.notebook("./kernel.jl", pwd()*"/../examples/"; execute)
Literate.markdown("./kernel.jl", pwd()*"/../docs/src/examples/kernel/"; execute)
