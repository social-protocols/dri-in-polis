instantiate:
   julia --project -e 'using Pkg; Pkg.instantiate()'

run-reference:
   julia --project reference_dri.jl

download-polis-data:
	git submodule update --init --recursive

run-polis CASE="vtaiwan.uberx" METHOD="" THRESHOLD="":
   just download-polis-data
   julia --project polis_dri.jl {{CASE}} {{METHOD}} {{THRESHOLD}}


run-random mode="delta" N="1000": 
   julia --project random_dri.jl  {{mode}} {{N}}