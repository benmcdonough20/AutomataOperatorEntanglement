#/bin/bash

for f in src/produce_figures/*.jl; do julia "$f"; done