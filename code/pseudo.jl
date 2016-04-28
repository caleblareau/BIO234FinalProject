## MCT compression pseudocode

import data

foreach i in ethnicity
    MCT[1,i] = mean(haplotype[1,ethnicity==i])
    foreach j in markers
        MCT[j,i] = mean(haplotype[j]==haplotype[j-1])
    end
end

export MCT, ethnicity_key

## MCT expansion pseudocode
import MCT, ethnicity_key

n=number_individuals
pop = matrix(NA,markers,n)
eth = random(ethnicity,n)
foreach i in ethnicity
    rand = random(0,1,length=n[eth==i])
    start.state = rand <= MCT[1,i[eth==1]]
    foreach j in markers
        rand = random(0,1,length=n[eth==i])
        transitions[j] = rand <= MCT[j,i[eth==1]]
    end
end

export transitions, ethnicity_key


