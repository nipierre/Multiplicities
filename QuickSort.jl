using DelimitedFiles

HadronN = readdlm("test.txt")
HadronM = readdlm("out1.txt")
list = Int64[]

for i in 1:4280
    if (((HadronN[i,7] == 0 || HadronN[i,7] == 1) && HadronM[i,7] == 0)
        || ((HadronN[i,7] == 2 || HadronN[i,7] == 3) && HadronM[i,7] == 1)
        || ((HadronN[i,7] == 4 || HadronN[i,7] == 5) && HadronM[i,7] == 2)
        || ((HadronN[i,7] == 6 || HadronN[i,7] == 7) && HadronM[i,7] == -1))
    else
        push!(list,i)
    end
end

writedlm("plots/QuickSort.txt", list)
