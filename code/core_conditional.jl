# get transition probabilities
function getTransitions(loc)
    idx = haplogroup.==loc;
    init = 1-sum(top[idx'])/sum(idx); #proportion of 0s in first position
    tP = 1-sum(sameState[:,idx[1:end]],2)./sum(idx); #probability of transitioning to a different state
    tP = tP[1:end-1]; #remove last row
    cat(1, init, tP); # add initial probability to top with transitions following
end

# Generic column means function
function colMeans(arr)
    out = mean(arr,1)
    return(out)
end

# make supergroup populations
function superGroup(loc)
    idx = haplogroup.==loc;
    if loc=="ASW" || loc=="LWK" || loc=="YRI"
        haplogroup[idx] = "AFR"
    elseif loc=="CEU" || loc=="FIN" || loc=="GBR" || loc == "IBS" || loc=="TSI"
        haplogroup[idx] = "EUR"
    elseif loc=="PUR" || loc=="MXL" || loc=="CLM"
        haplogroup[idx] = "AMR"
    elseif loc=="CHB" || loc=="CHS" || loc=="JPT"
        haplogroup[idx] = "EAS"
    end
end

# Thin samples
function thinTransitions(transMat, n)
    head = transMat[1,:];
    tmean = map(colMeans,[transMat[i:i+n-1,:] for i=2:n:size(transMat)[1]-n+1])
    out = vcat(tmean...)
    out = cat(1,head,out)
    writedlm(string(path,"output/transitionMat_t", n, ".txt"), out);
    return out;
end

#simulate population from transitions matrix
function simulatePopulation(transMat, levels, n, thin)
    sampinfo = levels[rand(1:end,n)]; # generate random ethnicity
    mat = Array(Int8,size(transMat)[1],2*n); #allocate empty array
    sampexpand = repeat(sampinfo,inner=[2])
    for k=1:length(levels)
        loc=levels[k]; 
        idx = sampexpand.==loc;
        # populate first row using population proportion of zeros
        r = rand(sum(idx));
        idz = 1-(r.<=transMat[1,k]);
        mat[1,idx]=idz;
        # populate subsequent rows using transition probabilities
        for i=2:size(mat)[1]
            r = rand(sum(idx));
            idz = r.<=transMat[i,k];
            for j=1:length(idz)
                if idz[j] == true
                    mat[i,find(idx)[j]] = 1-mat[i-1,find(idx)[j]]; #transitions states
                else 
                    mat[i,find(idx)[j]] = mat[i-1,find(idx)[j]]; #does not transition states
                end
            end
        end
    end
    writedlm(string(path,"input/sim_", n, "_t", thin,".hap"), mat);
    writedlm(string(path,"input/sim_", n, "_t", thin,".sample"), sampinfo);
end

# Create MCT matrix
# import data
hapdat = readdlm(string(path, hap), header=false);
sampdat = readdlm(string(path, samp), Any, skipstart=1);

# Get groups per haplotype; two haplotypes (columns) per sample
levels = sort(unique(sampdat[:,2]));
haplogroup = repeat(sampdat[:,2],inner=[2]);

# Compute Transition probabilities
top = copy(hapdat[1,:]);
hapdat_shift =  cat(1,copy(hapdat[2:end,:]),top); #shift first row to bottom
sameState = hapdat.==hapdat_shift;
sameState = convert(Array{Uint8}, bitpack(sameState));

map(superGroup,levels);
levels = sort(unique(haplogroup));
trans = map(getTransitions,levels);
trans = hcat(trans...);

writedlm(string(path,"output/transitionMat.txt"), trans);
writedlm(string(path,"output/legend.txt"), levels');
