instantiate:
   julia --project -e 'using Pkg; Pkg.instantiate()'

reference-implementation:
   julia --project reference-implementation.jl

download-polis-data:
   git submodule update --init --recursive

poc-polis METHOD="" THRESHOLD="":
   just download-polis-data
   julia --project poc-polis.jl {{METHOD}} {{THRESHOLD}}

random-tagging N="" mode="" : 
   julia --project random-tagging.jl {{N}} {{mode}}

validity-test:
   julia --project validity-test.jl

statement-subset: 
   julia --project statement-subset.jl
   