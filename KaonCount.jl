using DelimitedFiles

KaonN = readdlm("data/multiplicities_kaon_norich.txt")
KaonM = readdlm("data/Kaon_RC_VM_20190327_noRICHunfold.txt")

Count = zeros(540,9)

for i in 1:540
    Count[i,1] = KaonN[i,1]
    Count[i,2] = KaonN[i,2]
    Count[i,3] = KaonN[i,3]
    Count[i,4] = KaonN[i,11]-KaonM[i,10]
    Count[i,5] = KaonN[i,12]-KaonM[i,11]
    Count[i,6] = KaonN[i,21]-KaonM[i,12]
    Count[i,7] = KaonN[i,22]-KaonM[i,13]
    Count[i,8] = KaonN[i,9]>0 ? (KaonN[i,8]-KaonM[i,4])/KaonN[i,9] : 0
    Count[i,9] = KaonN[i,18]>0 ? (KaonN[i,18]-KaonM[i,6])/KaonN[i,19] : 0
end

writedlm("plots/CountKaon.txt", Count)
