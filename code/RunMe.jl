# path declarations
path = "/Users/Ploenzke/Documents/Harvard/Data\ Structures/Final/BIO234FinalProject/";
hap = "input/ALL_1000G_phase1interim_jun2011_chr20_impute.hap";
samp = "input/ALL_1000G_phase1interim_jun2011.sample";
legend = "input/ALL_1000G_phase1interim_jun2011_chr20_impute.legend";

# time compression
@time include(string(path,"code/core.jl"));

# thin samples
trans2 = thinTransitions(trans,2);
trans5 = thinTransitions(trans,5);
trans10 = thinTransitions(trans,10);
trans100 = thinTransitions(trans,100);

sampsize = 1094
@time simulatePopulation(trans, levels, sampsize, 0)
@time simulatePopulation(trans2, levels, sampsize, 2)
@time simulatePopulation(trans5, levels, sampsize, 5)
@time simulatePopulation(trans10, levels, sampsize, 10)
@time simulatePopulation(trans100, levels, sampsize, 100)
