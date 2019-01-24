using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]

DVMp = readdlm("data/DVM.dat")
DVMd = readdlm("data/DVM_2006.dat")
DVMa = readdlm("data/DVM_Andrea.dat")

DVMpipp = zeros((9,6,12))
DVMpimp = zeros((9,6,12))
DVMKpp = zeros((9,6,12))
DVMKmp = zeros((9,6,12))
DVMpipd = zeros((9,6,12))
DVMKpd = zeros((9,6,12))
DVMpimd = zeros((9,6,12))
DVMKmd = zeros((9,6,12))
DVMpipa = zeros((9,6,12))
DVMKpa = zeros((9,6,12))
DVMpima = zeros((9,6,12))
DVMKma = zeros((9,6,12))

for i in 1:648
    # DVMpipp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = DVMp[i,4]/DVMp[i,5]
    # DVMKpp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = DVMp[i,8]/DVMp[i,9]
    # DVMpimp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = DVMp[i,6]/DVMp[i,7]
    # DVMKmp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = DVMp[i,10]/DVMp[i,11]
    DVMpipp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = 1/(-DVMp[i,4]+1)-1
    DVMKpp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = 1/(-DVMp[i,8]+1)-1
    DVMpimp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = 1/(-DVMp[i,6]+1)-1
    DVMKmp[Int.(DVMp[i,1]),Int.(DVMp[i,2]),Int.(DVMp[i,3])] = 1/(-DVMp[i,10]+1)-1
end

for i in 1:540
    # DVMpipa[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = (1-DVMa[i,5])/(1-DVMa[i,4])
    # DVMKpa[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = (1-DVMa[i,7])/(1-DVMa[i,4])
    # DVMpima[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = (1-DVMa[i,6])/(1-DVMa[i,4])
    # DVMKma[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = (1-DVMa[i,8])/(1-DVMa[i,4])
    DVMpipa[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = 1/DVMa[i,5]-1
    DVMKpa[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = 1/DVMa[i,7]-1
    DVMpima[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = 1/DVMa[i,6]-1
    DVMKma[Int.(DVMa[i,1]+1),Int.(DVMa[i,2]+1),Int.(DVMa[i,3]+1)] = 1/DVMa[i,8]-1
end

for i in 1:310
    DVMpipd[Int.(DVMd[i,1]),Int.(DVMd[i,2]),Int.(DVMd[i,3])] = DVMd[i,4]/DVMd[i,5]
    DVMKpd[Int.(DVMd[i,1]),Int.(DVMd[i,2]),Int.(DVMd[i,3])] = DVMd[i,8]/DVMd[i,9]
    DVMpimd[Int.(DVMd[i,1]),Int.(DVMd[i,2]),Int.(DVMd[i,3])] = DVMd[i,6]/DVMd[i,7]
    DVMKmd[Int.(DVMd[i,1]),Int.(DVMd[i,2]),Int.(DVMd[i,3])] = DVMd[i,10]/DVMd[i,11]
end

zmid = zeros(12)
for i in 1:12
    zmid[i] = (z[i+1]+z[i])/2
end

DVMpip = zeros((9,6,12))
DVMpim = zeros((9,6,12))
DVMKp = zeros((9,6,12))
DVMKm = zeros((9,6,12))
DVMpipred = zeros(12,45)
DVMKpred = zeros(12,45)
DVMpimred = zeros(12,45)
DVMKmred = zeros(12,45)

for i in 1:9
    for j in 1:5
        global DVMpip, DVMpim, DVMKp, DVMKm
        global DVMpired, DVMKred
        for k in 1:12
            DVMpipred[k,i+(j-1)*9] = (DVMpipp[i,j,k] > 0 && DVMpipa[i,j,k] > 0) ? DVMpipp[i,j,k]/DVMpipa[i,j,k] : 0
            DVMpimred[k,i+(j-1)*9] = (DVMpimp[i,j,k] > 0 && DVMpima[i,j,k] > 0) ? DVMpimp[i,j,k]/DVMpima[i,j,k] : 0
            DVMKpred[k,i+(j-1)*9] = (DVMKpp[i,j,k] > 0 && DVMKpa[i,j,k] > 0) ? DVMKpp[i,j,k]/DVMKpa[i,j,k] : 0
            DVMKmred[k,i+(j-1)*9] = (DVMKmp[i,j,k] > 0 && DVMKma[i,j,k] > 0) ? DVMKmp[i,j,k]/DVMKma[i,j,k] : 0
            # DVMpipred[k,i+(j-1)*9] = (DVMpipp[i,j,k] > 0 && DVMpipd[i,j,k] > 0) ? DVMpipp[i,j,k]/DVMpipd[i,j,k] : 0
            # DVMpimred[k,i+(j-1)*9] = (DVMpimp[i,j,k] > 0 && DVMpimd[i,j,k] > 0) ? DVMpimp[i,j,k]/DVMpimd[i,j,k] : 0
            # DVMKpred[k,i+(j-1)*9] = (DVMKpp[i,j,k] > 0 && DVMKpd[i,j,k] > 0) ? DVMKpp[i,j,k]/DVMKpd[i,j,k] : 0
            # DVMKmred[k,i+(j-1)*9] = (DVMKmp[i,j,k] > 0 && DVMKmd[i,j,k] > 0) ? DVMKmp[i,j,k]/DVMKmd[i,j,k] : 0
        end
    end
end

plot(zmid,DVMpipred,layout=(5,9),seriestype=:scatter,
                                 size = (3200,1800),
                                 dpi = 300,
                                 ylims = (0.,2.),
                                 # ylims = (0.75,1.25),
                                 legend=false)
plot!(zmid,DVMpimred,layout=(5,9),seriestype=:scatter,
                                  size = (3200,1800),
                                  dpi = 300,
                                  ylims = (0.,2.),
                                  # ylims = (0.75,1.25),
                                  legend=false)
savefig("plots/DVMComppi.png")

plot(zmid,DVMKpred,layout=(5,9),seriestype=:scatter,
                                size = (3200,1800),
                                dpi = 300,
                                ylims = (0.,2.),
                                # ylims = (0.75,1.25),
                                legend=false)
plot!(zmid,DVMKmred,layout=(5,9),seriestype=:scatter,
                                 size = (3200,1800),
                                 dpi = 300,
                                 ylims = (0.,2.),
                                 # ylims = (0.75,1.25),
                                 legend=false)
savefig("plots/DVMCompK.png")
