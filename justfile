instantiate:
   julia --project -e 'using Pkg; Pkg.instantiate()'

reference-implementation:
   julia --project reference-implementation.jl

download-polis-data:
   git submodule update --init --recursive

poc-polis CASE="vtaiwan.uberx" METHOD="" THRESHOLD="25":
   just download-polis-data
   julia --project poc-polis.jl {{CASE}} {{METHOD}} {{THRESHOLD}}

random-tagging N="" mode="" : 
   julia --project random-tagging.jl {{N}} {{mode}}

statement-subset: 
   julia --project statement-subset.jl
   