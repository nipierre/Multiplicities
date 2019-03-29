using DelimitedFiles

HadronN = readdlm("test.txt")
HadronM = readdlm("out1.txt")
list = Int64[]

for i in 1:4280
    if (abs(HadronN[i,1]-HadronM[i,1])> 0.00001 || abs(HadronN[i,2]-HadronM[i,2]) > 0.01 || abs(HadronN[i,3]-HadronM[i,3]) > 0.0001
        || abs(HadronN[i,4]-HadronM[i,4]) > 0.01 || abs(HadronN[i,5]-HadronM[i,5]) > 0.0001)
        push!(list,i)
    end
end

writedlm("plots/QuickSort.txt", list)
