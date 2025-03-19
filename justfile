instantiate:
   julia --project -e 'using Pkg; Pkg.instantiate()'

run-reference:
   julia --project reference_dri.jl

download-polis-data:
	git submodule update --init --recursive

run-polis CASE="vtaiwan.uberx" METHOD="phi":
   just download-polis-data
   julia --project polis_dri.jl {{CASE}} {{METHOD}}
