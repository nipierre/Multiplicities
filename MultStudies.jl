using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]
xexp = [.007,.015,.025,.035,.05,.08,.12,.16,.29]
zFFred = [3,8,13,18,23,28,33,38,43,48,53,61]

Mult1 = readdlm("data/Mult2016Hadron.txt")
Mult2 = readdlm("data/Mult2016HadronTheta.txt")

Multp1 = zeros((9,12))
Multm1 = zeros((9,12))
Multp2 = zeros((9,12))
Multm2 = zeros((9,12))
sigp1 = zeros((9,12))
sigm1 = zeros((9,12))
sigp2 = zeros((9,12))
sigm2 = zeros((9,12))


for i in 1:9
    for j in 1:12
        global Multm1, Multm2, sigm1 , sigm2
        Multm1[i,j] = Mult1[j+24*(i-1),4]
        Multm2[i,j] = Mult2[j+24*(i-1),4]
        sigm1[i,j] = Mult1[j+24*(i-1),5]
        sigm2[i,j] = Mult2[j+24*(i-1),5]
    end
    for j in 1:12
        global Multp1, Multp2, sigp1 , sigp2
        Multp1[i,j] = Mult1[j+24*(i-1)+12,4]
        Multp2[i,j] = Mult2[j+24*(i-1)+12,4]
        sigp1[i,j] = Mult1[j+24*(i-1),5]
        sigp2[i,j] = Mult2[j+24*(i-1),5]
    end
end

zmid = zeros(12)
for i in 1:12
    zmid[i] = (z[i+1]+z[i])/2
end

Ms1 = zeros(9)
Ms2 = zeros(9)
Ss1 = zeros(9)
Ss2 = zeros(9)
Mp1 = zeros(12)
Mp2 = zeros(12)
Sp1 = zeros(12)
Sp2 = zeros(12)
Mr1 = zeros(9)
Mr2 = zeros(9)
Sr1 = zeros(9)
Sr2 = zeros(9)
Mm1 = zeros(12)
Mm2 = zeros(12)
Sm1 = zeros(12)
Sm2 = zeros(12)
Mrp1 = zeros(9)
Mrp2 = zeros(9)
Mrm1 = zeros(9)
Mrm2 = zeros(9)
Srp1 = zeros(9)
Srp2 = zeros(9)
Srm1 = zeros(9)
Srm2 = zeros(9)

flag = 0

for i in 1:9
    global Mr1, Mr2
    global Sr1, Sr2
    global Ms1, Ms2, Mp1, Mp2, Mm1, Mm2, Mrp1, Mrp2, Mrm1, Mrm2
    global Ss1, Ss2, Sp1, Sp2, Sm1, Sm2, Srp1, Srp2, Srm1, Srm2
    global flag
    flag = 0
    for j in 1:12
        if flag < 1
            flag = (Multp2[i,j] > 0 && Multm2[i,j] > 0) ? 0 : 1
        end
        Ms1[i] += (Multp1[i,j]+Multm1[i,j])*(z[j+1]-z[j])
        Ms2[i] += (Multp2[i,j]+Multm2[i,j])*(z[j+1]-z[j])
        Ss1[i] += (sigp1[i,j]+sigm1[i,j])*(z[j+1]-z[j])
        Ss2[i] += (sigp2[i,j]+sigm2[i,j])*(z[j+1]-z[j])
        Mrp1[i] += Multp1[i,j]*(z[j+1]-z[j])
        Mrp2[i] += Multp2[i,j]*(z[j+1]-z[j])
        Mrm1[i] += Multm1[i,j]*(z[j+1]-z[j])
        Mrm2[i] += Multm2[i,j]*(z[j+1]-z[j])
        Srp1[i] += sigp1[i,j]*(z[j+1]-z[j])
        Srp2[i] += sigp2[i,j]*(z[j+1]-z[j])
        Srm1[i] += sigm1[i,j]*(z[j+1]-z[j])
        Srm2[i] += sigm2[i,j]*(z[j+1]-z[j])
        Mp1[j] = Multp1[i,j]
        Mm1[j] = Multm1[i,j]
        Mp2[j] = Multp2[i,j]
        Mm2[j] = Multm2[i,j]
        Sp1[j] = sqrt(sigp1[i,j])
        Sm1[j] = sqrt(sigm1[i,j])
        Sp2[j] = sqrt(sigp2[i,j])
        Sm2[j] = sqrt(sigm2[i,j])
    end
    t = string("Pion Multiplicities at x =", xexp[i])
    plot(zmid, Mp1, lw=3, xlims = (0.,1),
                           ylims = (-0.1,3),
                           linecolor = :red,
                           xlabel = "z",
                           ylabel = L"M^{h}",
                           ribbon = Sp1,
                           label = L"h^+_{RD}")
    plot!(zmid, Mm1, lw=3, xlims = (0.,1),
                            ylims = (-0.1,3),
                            linecolor = :blue,
                            xlabel = "z",
                            ylabel = L"M^{h}",
                            ribbon = Sm1,
                            label = L"h^-_{RD}")
    plot!(zmid, Mp2, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :red,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{h}",
                             ribbon = Sp2,
                             label = L"h^+_{MC}")
    plot!(zmid, Mm2, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :blue,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{h}",
                             ribbon = Sm2,
                             label = L"h^-_{MC}",
                             title=t)
    savefig(string("plots/MultiplicitiesTest",i,".png"))
    Mr1[i] = (Mrp1[i]/Mrm1[i])
    Mr2[i] = (Mrp2[i]/Mrm2[i])
    Mr1[i] = (Mrp1[i]/Mrm1[i])
    Mr2[i] = (Mrp2[i]/Mrm2[i])
    Sr1[i] = sqrt(Srp1[i]+Srm1[i])
    Sr2[i] = sqrt(Srp2[i]+Srm2[i])
    Sr1[i] = sqrt(Srp1[i]+Srm1[i])
    Sr2[i] = sqrt(Srp2[i]+Srm2[i])
    Ss1[i] = sqrt(Ss1[i])
    Ss2[i] = sqrt(Ss2[i])
    if flag > 0
        Ms2[i] = 0
        Ss2[i] = 0
        Mr2[i] = 0
        Sr2[i] = 0
    end
end

scatter(xexp,Ms1, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           ylims = (0.6,1.2),
           xlabel = "x",
           ylabel = L"\int M^{h^+}+M^{h^-} dz",
           yerror = Ss1,
           label = "RD")
scatter!(xexp,Ms2, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.6,1.2),
            xlabel = "x",
            ylabel = L"\int M^{h^+}+M^{h^-} dz",
            yerror = Ss2,
            label = "MC")
savefig("plots/MultiplicitiesSumTest.png")

scatter(xexp,Mr1, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           ylims = (1.1,1.9),
           xlabel = "x",
           ylabel = L"\frac{\int M^{h^+} dz}{\int M^{h^-} dz}",
           yerror = Sr1,
           label = "RD")
scatter!(xexp,Mr2, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (1.1,1.9),
            xlabel = "x",
            ylabel = L"\frac{\int M^{h^+} dz}{\int M^{h^-} dz}",
            yerror = Sr2,
            label = "MC")
savefig("plots/MultiplicitiesRatioTest.png")
