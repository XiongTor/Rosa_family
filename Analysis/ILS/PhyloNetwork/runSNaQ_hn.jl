#!/usr/bin/env julia

# file "runSNaQ.jl". run in the shell like this in general:
# julia runSNaQ.jl hvalue nruns
# example for h=2 and default 10 runs:
# julia runSNaQ.jl 2
# or example for h=3 and 50 runs:
# julia runSNaQ.jl 3 50

length(ARGS) > 0 ||
    error("need 1 or 2 arguments: # reticulations (h) and # runs (optional, 10 by default)")
h = parse(Int, ARGS[1])
nruns = 10
if length(ARGS) > 1
    nruns = parse(Int, ARGS[2])
end
outputfile = string("net", h, "_", nruns, "runs") # example: "net2_10runs"
seed = 1234 + h # change as desired! Best to have it different for different h
@info "will run SNaQ with h=$h, # of runs=$nruns, seed=$seed, output will go to: $outputfile"

using Distributed
addprocs(nruns)
@everywhere using PhyloNetworks
net0_h6 = readTopology("astral.tre");  #读取起始树，为了避免并行时linux系统环境变量得区分，在h为1时设置为net0_h1
using DataFrames, CSV
df_sp = DataFrame(CSV.File("tableCF.csv", pool=false); copycols=false); #读取CF表
d_sp = readTableCF!(df_sp);
net_h1 = snaq!(net0_h6, d_sp, hmax=h, filename=outputfile, seed=seed, runs=nruns) 