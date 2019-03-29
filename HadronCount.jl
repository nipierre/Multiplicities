using DelimitedFiles

HadronN = readdlm("data/multiplicities_kaon_norich.txt")
HadronM = readdlm("data/Kaon_RC_VM_20190327_noRICHunfold.txt")

PionN = readdlm("data/multiplicities_kaon_norich.txt")
PionM = readdlm("data/Kaon_RC_VM_20190327_noRICHunfold.txt")

KaonN = readdlm("data/multiplicities_kaon_norich.txt")
KaonM = readdlm("data/Kaon_RC_VM_20190327_noRICHunfold.txt")

ProtonN = readdlm("data/multiplicities_kaon_norich.txt")
ProtonM = readdlm("data/Kaon_RC_VM_20190327_noRICHunfold.txt")

CountH = zeros(540,9)
CountP = zeros(540,9)
CountK = zeros(540,9)
CountPr = zeros(540,9)

for i in 1:540
    CountH[i,1] = HadronN[i,1]
    CountH[i,2] = HadronN[i,2]
    CountH[i,3] = HadronN[i,3]
    CountH[i,4] = HadronN[i,11]-HadronM[i,10]
    CountH[i,5] = HadronN[i,12]-HadronM[i,11]
    CountH[i,6] = HadronN[i,21]-HadronM[i,12]
    CountH[i,7] = HadronN[i,22]-HadronM[i,13]
    CountH[i,8] = HadronN[i,9]>0 ? (HadronN[i,8]-HadronM[i,4])/HadronN[i,9] : 0
    CountH[i,9] = HadronN[i,18]>0 ? (HadronN[i,18]-HadronM[i,6])/HadronN[i,19] : 0

    CountP[i,1] = PionN[i,1]
    CountP[i,2] = PionN[i,2]
    CountP[i,3] = PionN[i,3]
    CountP[i,4] = PionN[i,11]-PionM[i,10]
    CountP[i,5] = PionN[i,12]-PionM[i,11]
    CountP[i,6] = PionN[i,21]-PionM[i,12]
    CountP[i,7] = PionN[i,22]-PionM[i,13]
    CountP[i,8] = PionN[i,9]>0 ? (PionN[i,8]-PionM[i,4])/PionN[i,9] : 0
    CountP[i,9] = PionN[i,18]>0 ? (PionN[i,18]-PionM[i,6])/PionN[i,19] : 0

    CountK[i,1] = KaonN[i,1]
    CountK[i,2] = KaonN[i,2]
    CountK[i,3] = KaonN[i,3]
    CountK[i,4] = KaonN[i,11]-KaonM[i,10]
    CountK[i,5] = KaonN[i,12]-KaonM[i,11]
    CountK[i,6] = KaonN[i,21]-KaonM[i,12]
    CountK[i,7] = KaonN[i,22]-KaonM[i,13]
    CountK[i,8] = KaonN[i,9]>0 ? (KaonN[i,8]-KaonM[i,4])/KaonN[i,9] : 0
    CountK[i,9] = KaonN[i,18]>0 ? (KaonN[i,18]-KaonM[i,6])/KaonN[i,19] : 0

    CountPr[i,1] = ProtonN[i,1]
    CountPr[i,2] = ProtonN[i,2]
    CountPr[i,3] = ProtonN[i,3]
    CountPr[i,4] = ProtonN[i,11]-ProtonM[i,10]
    CountPr[i,5] = ProtonN[i,12]-ProtonM[i,11]
    CountPr[i,6] = ProtonN[i,21]-ProtonM[i,12]
    CountPr[i,7] = ProtonN[i,22]-ProtonM[i,13]
    CountPr[i,8] = ProtonN[i,9]>0 ? (ProtonN[i,8]-ProtonM[i,4])/ProtonN[i,9] : 0
    CountPr[i,9] = ProtonN[i,18]>0 ? (ProtonN[i,18]-ProtonM[i,6])/ProtonN[i,19] : 0
end

writedlm("plots/CountHadron.txt", CountH)
writedlm("plots/CountPion.txt", CountP)
writedlm("plots/CountKaon.txt", CountK)
writedlm("plots/CountProton.txt", CountPr)
