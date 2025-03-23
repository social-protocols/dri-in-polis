instantiate:
   julia --project -e 'using Pkg; Pkg.instantiate()'

run-reference:
   julia --project reference_dri.jl

download-polis-data:
	git submodule update --init --recursive

run-polis CASE="vtaiwan.uberx" METHOD="" THRESHOLD="":
   just download-polis-data
   julia --project polis_dri.jl {{CASE}} {{METHOD}} {{THRESHOLD}}


run-resampling N="" mode="" : 
   julia --project resampling_dri.jl {{N}} {{mode}}