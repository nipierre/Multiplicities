using DelimitedFiles

Hadronp = readdlm("data/h+deut.txt")
DVMhp = readdlm("data/h+DVM.txt")
RChp = readdlm("data/h+RC.txt")
Hadronm = readdlm("data/h-deut.txt")
DVMhm = readdlm("data/h-DVM.txt")
RChm = readdlm("data/h-RC.txt")
Pionp = readdlm("data/pi+deut.txt")
DVMpp = readdlm("data/pi+DVM.txt")
RCpp = readdlm("data/pi+RC.txt")
Pionm = readdlm("data/pi-deut.txt")
DVMpm = readdlm("data/pi-DVM.txt")
RCpm = readdlm("data/pi-RC.txt")
DVMDIS = readdlm("data/DISDVM.txt")
RCDIS = readdlm("data/DISRC.txt")

Mh = zeros(311,19)
Mpi = zeros(311,19)

for i in 1:311
    Mh[i,1] = Hadronp[i,2]
    Mh[i,2] = Hadronp[i,5]
    Mh[i,3] = Hadronp[i,9]
    Mh[i,4] = Hadronp[i,1]
    Mh[i,5] = Hadronp[i,4]
    Mh[i,6] = Hadronp[i,7]
    Mh[i,7] = Hadronp[i,8]
    Mh[i,8] = Hadronp[i,11]*(DVMhp[i,11]/DVMDIS[i,11])*(RChp[i,11]/RCDIS[i,11])
    Mh[i,9] = Hadronp[i,12]
    Mh[i,10] = Hadronp[i,14]
    Mh[i,11] = Mh[i,8]>0 ? 1 : 0
    Mh[i,12] = Hadronm[i,1]
    Mh[i,13] = Hadronm[i,4]
    Mh[i,14] = Hadronm[i,7]
    Mh[i,15] = Hadronm[i,8]
    Mh[i,16] = Hadronm[i,11]*(DVMhm[i,11]/DVMDIS[i,11])*(RChm[i,11]/RCDIS[i,11])
    Mh[i,17] = Hadronm[i,12]
    Mh[i,18] = Hadronm[i,14]
    Mh[i,19] = Mh[i,16]>0 ? 1 : 0

    Mpi[i,1] = Pionp[i,2]
    Mpi[i,2] = Pionp[i,5]
    Mpi[i,3] = Pionp[i,9]
    Mpi[i,4] = Pionp[i,1]
    Mpi[i,5] = Pionp[i,4]
    Mpi[i,6] = Pionp[i,7]
    Mpi[i,7] = Pionp[i,8]
    Mpi[i,8] = Pionp[i,11]*(DVMpp[i,11]/DVMDIS[i,11])*(RCpp[i,11]/RCDIS[i,11])
    Mpi[i,9] = Pionp[i,12]
    Mpi[i,10] = Pionp[i,14]
    Mpi[i,11] = Mpi[i,8]>0 ? 1 : 0
    Mpi[i,12] = Pionm[i,1]
    Mpi[i,13] = Pionm[i,4]
    Mpi[i,14] = Pionm[i,7]
    Mpi[i,15] = Pionm[i,8]
    Mpi[i,16] = Pionm[i,11]*(DVMpm[i,11]/DVMDIS[i,11])*(RCpm[i,11]/RCDIS[i,11])
    Mpi[i,17] = Pionm[i,12]
    Mpi[i,18] = Pionm[i,14]
    Mpi[i,19] = Mpi[i,16]>0 ? 1 : 0
end

writedlm("data/MultiplicityHadron_2006.txt", Mh)
writedlm("data/MultiplicityPion_2006.txt", Mpi)
