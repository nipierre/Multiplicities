using DelimitedFiles
using Plots
using LaTeXStrings

z = [.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85]
xexp = [.007,.015,.025,.035,.05,.08,.12,.16,.29]
zFFred = [3,8,13,18,23,28,33,38,43,48,53,61]

Mult = readdlm("data/Mult2016HadronVertexedMU+.txt")

Multp1 = zeros((9,12))
Multm1 = zeros((9,12))
Multp2 = zeros((9,12))
Multm2 = zeros((9,12))
Multp3 = zeros((9,12))
Multm3 = zeros((9,12))
Multp4 = zeros((9,12))
Multm4 = zeros((9,12))
sigp1 = zeros((9,12))
sigm1 = zeros((9,12))
sigp2 = zeros((9,12))
sigm2 = zeros((9,12))
sigp3 = zeros((9,12))
sigm3 = zeros((9,12))
sigp4 = zeros((9,12))
sigm4 = zeros((9,12))


for i in 1:9
    for j in 1:12
        global Multm1, Multm2, sigm1 , sigm2
        Multm1[i,j] = Mult[j+24*(i-1),4]
        Multm2[i,j] = Mult[j+24*(i-1),6]
        Multm3[i,j] = Mult[j+24*(i-1),8]
        Multm4[i,j] = Mult[j+24*(i-1),10]
        sigm1[i,j] = Mult[j+24*(i-1),5]
        sigm2[i,j] = Mult[j+24*(i-1),7]
        sigm3[i,j] = Mult[j+24*(i-1),9]
        sigm4[i,j] = Mult[j+24*(i-1),11]
    end
    for j in 1:12
        global Multp1, Multp2, sigp1 , sigp2
        Multp1[i,j] = Mult[j+24*(i-1)+12,4]
        Multp2[i,j] = Mult[j+24*(i-1)+12,6]
        Multp3[i,j] = Mult[j+24*(i-1)+12,8]
        Multp4[i,j] = Mult[j+24*(i-1)+12,10]
        sigp1[i,j] = Mult[j+24*(i-1),5]
        sigp2[i,j] = Mult[j+24*(i-1),7]
        sigp3[i,j] = Mult[j+24*(i-1),9]
        sigp4[i,j] = Mult[j+24*(i-1),11]
    end
end

zmid = zeros(12)
for i in 1:12
    zmid[i] = (z[i+1]+z[i])/2
end

Ms1 = zeros(9)
Ms2 = zeros(9)
Ms3 = zeros(9)
Ms4 = zeros(9)
Ss1 = zeros(9)
Ss2 = zeros(9)
Ss3 = zeros(9)
Ss4 = zeros(9)
Mp1 = zeros(12)
Mp2 = zeros(12)
Mp3 = zeros(12)
Mp4 = zeros(12)
Sp1 = zeros(12)
Sp2 = zeros(12)
Sp3 = zeros(12)
Sp4 = zeros(12)
Mr1 = zeros(9)
Mr2 = zeros(9)
Mr3 = zeros(9)
Mr4 = zeros(9)
Sr1 = zeros(9)
Sr2 = zeros(9)
Sr3 = zeros(9)
Sr4 = zeros(9)
Mm1 = zeros(12)
Mm2 = zeros(12)
Mm3 = zeros(12)
Mm4 = zeros(12)
Sm1 = zeros(12)
Sm2 = zeros(12)
Sm3 = zeros(12)
Sm4 = zeros(12)
Mrp1 = zeros(9)
Mrp2 = zeros(9)
Mrp3 = zeros(9)
Mrp4 = zeros(9)
Mrm1 = zeros(9)
Mrm2 = zeros(9)
Mrm3 = zeros(9)
Mrm4 = zeros(9)
Srp1 = zeros(9)
Srp2 = zeros(9)
Srp3 = zeros(9)
Srp4 = zeros(9)
Srm1 = zeros(9)
Srm2 = zeros(9)
Srm3 = zeros(9)
Srm4 = zeros(9)

flag1 = 0
flag2 = 0
flag3 = 0
flag4 = 0

for i in 1:9
    global Mr1, Mr2, Mr3, Mr4
    global Sr1, Sr2, Sr3, Sr4
    global Ms1, Ms2, Mp1, Mp2, Mm1, Mm2, Mrp1, Mrp2, Mrm1, Mrm2
    global Ms3, Ms4, Mp3, Mp4, Mm3, Mm4, Mrp3, Mrp4, Mrm3, Mrm4
    global Ss1, Ss2, Sp1, Sp2, Sm1, Sm2, Srp1, Srp2, Srm1, Srm2
    global Ss3, Ss4, Sp3, Sp4, Sm3, Sm4, Srp3, Srp4, Srm3, Srm4
    global flag1, flag2, flag3, flag4
    flag1 = 0
    flag2 = 0
    flag3 = 0
    flag4 = 0
    for j in 1:12
        if flag1 < 1
            flag1 = (Multp1[i,j] > 0 && Multm1[i,j] > 0) ? 0 : 1
        end
        if flag2 < 1
            flag2 = (Multp2[i,j] > 0 && Multm2[i,j] > 0) ? 0 : 1
        end
        if flag3 < 1
            flag3 = (Multp3[i,j] > 0 && Multm3[i,j] > 0) ? 0 : 1
        end
        if flag4 < 1
            flag4 = (Multp4[i,j] > 0 && Multm4[i,j] > 0) ? 0 : 1
        end
        Ms1[i] += (Multp1[i,j]+Multm1[i,j])*(z[j+1]-z[j])
        Ms2[i] += (Multp2[i,j]+Multm2[i,j])*(z[j+1]-z[j])
        Ms3[i] += (Multp3[i,j]+Multm3[i,j])*(z[j+1]-z[j])
        Ms4[i] += (Multp4[i,j]+Multm4[i,j])*(z[j+1]-z[j])
        Ss1[i] += (sigp1[i,j]+sigm1[i,j])*(z[j+1]-z[j])
        Ss2[i] += (sigp2[i,j]+sigm2[i,j])*(z[j+1]-z[j])
        Ss3[i] += (sigp3[i,j]+sigm3[i,j])*(z[j+1]-z[j])
        Ss4[i] += (sigp4[i,j]+sigm4[i,j])*(z[j+1]-z[j])
        Mrp1[i] += Multp1[i,j]*(z[j+1]-z[j])
        Mrp2[i] += Multp2[i,j]*(z[j+1]-z[j])
        Mrp3[i] += Multp3[i,j]*(z[j+1]-z[j])
        Mrp4[i] += Multp4[i,j]*(z[j+1]-z[j])
        Mrm1[i] += Multm1[i,j]*(z[j+1]-z[j])
        Mrm2[i] += Multm2[i,j]*(z[j+1]-z[j])
        Mrm3[i] += Multm3[i,j]*(z[j+1]-z[j])
        Mrm4[i] += Multm4[i,j]*(z[j+1]-z[j])
        Srp1[i] += sigp1[i,j]*(z[j+1]-z[j])
        Srp2[i] += sigp2[i,j]*(z[j+1]-z[j])
        Srp3[i] += sigp3[i,j]*(z[j+1]-z[j])
        Srp4[i] += sigp4[i,j]*(z[j+1]-z[j])
        Srm1[i] += sigm1[i,j]*(z[j+1]-z[j])
        Srm2[i] += sigm2[i,j]*(z[j+1]-z[j])
        Srm3[i] += sigm3[i,j]*(z[j+1]-z[j])
        Srm4[i] += sigm4[i,j]*(z[j+1]-z[j])
        Mp1[j] = Multp1[i,j]
        Mm1[j] = Multm1[i,j]
        Mp2[j] = Multp2[i,j]
        Mm2[j] = Multm2[i,j]
        Mp3[j] = Multp3[i,j]
        Mm3[j] = Multm3[i,j]
        Mp4[j] = Multp4[i,j]
        Mm4[j] = Multm4[i,j]
        Sp1[j] = sqrt(sigp1[i,j])
        Sm1[j] = sqrt(sigm1[i,j])
        Sp2[j] = sqrt(sigp2[i,j])
        Sm2[j] = sqrt(sigm2[i,j])
        Sp3[j] = sqrt(sigp3[i,j])
        Sm3[j] = sqrt(sigm3[i,j])
        Sp4[j] = sqrt(sigp4[i,j])
        Sm4[j] = sqrt(sigm4[i,j])
    end
    t = string("Pion Multiplicities at x =", xexp[i])
    plot(zmid, Mp1, lw=3, xlims = (0.,1),
                           ylims = (-0.1,3),
                           linecolor = :red,
                           xlabel = "z",
                           ylabel = L"M^{h}",
                           ribbon = Sp1,
                           label = L"h^+_{1}")
    plot!(zmid, Mm1, lw=3, xlims = (0.,1),
                            ylims = (-0.1,3),
                            linecolor = :blue,
                            xlabel = "z",
                            ylabel = L"M^{h}",
                            ribbon = Sm1,
                            label = L"h^-_{1}")
    plot!(zmid, Mp2, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :red,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{h}",
                             ribbon = Sp2,
                             label = L"h^+_{2}")
    plot!(zmid, Mm2, lw=3, xlims = (0.,1),
                             ylims = (-0.1,3),
                             linecolor = :blue,
                             linestyle = :dash,
                             xlabel = "z",
                             ylabel = L"M^{h}",
                             ribbon = Sm2,
                             label = L"h^-_{2}",
                             title=t)
     plot!(zmid, Mp3, lw=3, xlims = (0.,1),
                              ylims = (-0.1,3),
                              linecolor = :red,
                              linestyle = :dot,
                              xlabel = "z",
                              ylabel = L"M^{h}",
                              ribbon = Sp3,
                              label = L"h^+_{3}")
     plot!(zmid, Mm3, lw=3, xlims = (0.,1),
                              ylims = (-0.1,3),
                              linecolor = :blue,
                              linestyle = :dot,
                              xlabel = "z",
                              ylabel = L"M^{h}",
                              ribbon = Sm3,
                              label = L"h^-_{3}",
                              title=t)
      plot!(zmid, Mp4, lw=3, xlims = (0.,1),
                               ylims = (-0.1,3),
                               linecolor = :red,
                               linestyle = :dashdot,
                               xlabel = "z",
                               ylabel = L"M^{h}",
                               ribbon = Sp4,
                               label = L"h^+_{4}")
      plot!(zmid, Mm4, lw=3, xlims = (0.,1),
                               ylims = (-0.1,3),
                               linecolor = :blue,
                               linestyle = :dashdot,
                               xlabel = "z",
                               ylabel = L"M^{h}",
                               ribbon = Sm4,
                               label = L"h^-_{4}",
                               title=t)
    savefig(string("plots/MultiplicitiesVertex",i,".png"))
    Mr1[i] = (Mrp1[i]/Mrm1[i])
    Mr2[i] = (Mrp2[i]/Mrm2[i])
    Mr3[i] = (Mrp3[i]/Mrm3[i])
    Mr4[i] = (Mrp4[i]/Mrm4[i])
    Mr1[i] = (Mrp1[i]/Mrm1[i])
    Mr2[i] = (Mrp2[i]/Mrm2[i])
    Mr3[i] = (Mrp3[i]/Mrm3[i])
    Mr4[i] = (Mrp4[i]/Mrm4[i])
    Sr1[i] = sqrt(Srp1[i]+Srm1[i])
    Sr2[i] = sqrt(Srp2[i]+Srm2[i])
    Sr3[i] = sqrt(Srp3[i]+Srm3[i])
    Sr4[i] = sqrt(Srp4[i]+Srm4[i])
    Sr1[i] = sqrt(Srp1[i]+Srm1[i])
    Sr2[i] = sqrt(Srp2[i]+Srm2[i])
    Sr3[i] = sqrt(Srp3[i]+Srm3[i])
    Sr4[i] = sqrt(Srp4[i]+Srm4[i])
    Ss1[i] = sqrt(Ss1[i])
    Ss2[i] = sqrt(Ss2[i])
    Ss3[i] = sqrt(Ss3[i])
    Ss4[i] = sqrt(Ss4[i])
    if flag1 > 0
        Ms1[i] = 0
        Ss1[i] = 0
        Mr1[i] = 0
        Sr1[i] = 0
    end
    if flag2 > 0
        Ms2[i] = 0
        Ss2[i] = 0
        Mr2[i] = 0
        Sr2[i] = 0
    end
    if flag3 > 0
        Ms3[i] = 0
        Ss3[i] = 0
        Mr3[i] = 0
        Sr3[i] = 0
    end
    if flag4 > 0
        Ms4[i] = 0
        Ss4[i] = 0
        Mr4[i] = 0
        Sr4[i] = 0
    end
end

scatter(xexp,Ms1, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           ylims = (0.6,1.2),
           xlabel = "x",
           ylabel = L"\int M^{h^+}+M^{h^-} dz",
           yerror = Ss1,
           label = "1")
scatter!(xexp,Ms2, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.6,1.2),
            xlabel = "x",
            ylabel = L"\int M^{h^+}+M^{h^-} dz",
            yerror = Ss2,
            label = "2")
scatter!(xexp,Ms3, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.6,1.2),
            xlabel = "x",
            ylabel = L"\int M^{h^+}+M^{h^-} dz",
            yerror = Ss3,
            label = "3")
scatter!(xexp,Ms4, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (0.6,1.2),
            xlabel = "x",
            ylabel = L"\int M^{h^+}+M^{h^-} dz",
            yerror = Ss4,
            label = "4")
savefig("plots/MultiplicitiesSumVertex.png")

scatter(xexp,Mr1, lw=3,
           xscale = :log10,
           xlims = (0.01,1),
           ylims = (1.1,1.9),
           xlabel = "x",
           ylabel = L"\frac{\int M^{h^+} dz}{\int M^{h^-} dz}",
           yerror = Sr1,
           label = "1")
scatter!(xexp,Mr2, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (1.1,1.9),
            xlabel = "x",
            ylabel = L"\frac{\int M^{h^+} dz}{\int M^{h^-} dz}",
            yerror = Sr2,
            label = "2")
scatter!(xexp,Mr3, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (1.1,1.9),
            xlabel = "x",
            ylabel = L"\frac{\int M^{h^+} dz}{\int M^{h^-} dz}",
            yerror = Sr3,
            label = "3")
scatter!(xexp,Mr4, lw=3,
            xscale = :log10,
            xlims = (0.01,1),
            ylims = (1.1,1.9),
            xlabel = "x",
            ylabel = L"\frac{\int M^{h^+} dz}{\int M^{h^-} dz}",
            yerror = Sr4,
            label = "4")
savefig("plots/MultiplicitiesRatioVertex.png")
