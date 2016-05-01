# path declarations
path = "/Users/Ploenzke/Documents/Harvard/Data\ Structures/Final/BIO234FinalProject/";
hap = "input/ALL_1000G_phase1interim_jun2011_chr20_impute.hap";
samp = "input/ALL_1000G_phase1interim_jun2011.sample";
legend = "input/ALL_1000G_phase1interim_jun2011_chr20_impute.legend";

# time compression
@time include(string(path,"code/core.jl"));

# time thinning
@time trans2 = thinTransitions(trans,2);
@time trans5 = thinTransitions(trans,5);
@time trans10 = thinTransitions(trans,10);
@time trans100 = thinTransitions(trans,100);

# time data regeneration
sampsize = 100
@time trans = readdlm(string(path,"output/transitionMat.txt"), header=false);
@time simulatePopulation(trans, levels, sampsize, 0)
@time trans2 = readdlm(string(path,"output/transitionMat_t2.txt"), header=false);
@time simulatePopulation(trans2, levels, sampsize, 2)
@time trans5 = readdlm(string(path,"output/transitionMat_t5.txt"), header=false);
@time simulatePopulation(trans5, levels, sampsize, 5)
@time trans10 = readdlm(string(path,"output/transitionMat_t10.txt"), header=false);
@time simulatePopulation(trans10, levels, sampsize, 10)
@time trans100 = readdlm(string(path,"output/transitionMat_t100.txt"), header=false);
@time simulatePopulation(trans100, levels, sampsize, 100)
