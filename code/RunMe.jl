# path declarations
path = "/Users/Ploenzke/Documents/Harvard/Data\ Structures/Final/BIO234FinalProject-master/";
hap = "input/ALL_1000G_phase1interim_jun2011_chr20_impute.hap";
samp = "input/ALL_1000G_phase1interim_jun2011.sample";
legend = "input/ALL_1000G_phase1interim_jun2011_chr20_impute.legend";

# time compression
@time include(string(path,"code/core.jl"));

# thin samples
trans10 = thinTransitions(trans,2);
trans10 = thinTransitions(trans,10);
trans50 = thinTransitions(trans,50);
trans100 = thinTransitions(trans,100);

sampsize = 1000000
@time haps, sample = simulatePopulation(trans, levels, sampsize)
