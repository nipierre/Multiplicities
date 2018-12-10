using DelimitedFiles
using Plots
using LaTeXStrings

theta = [0.01,0.04,0.12]
mom06 = [9.,10.,11.,12.,13.,15.,17.,19.,22.,25.,27.,30.,35.,40.]
mom11 = [10.,11.,12.,13.,15.,17.,19.,22.,25.,27.,30.,35.,40.]
mom16 = [3.,5.,7.,10.,12.,13.,15.,17.,19.,22.,25.,27.,30.,35.,40.]

rich_2006_t1 = readdlm("data/rich_mat_2006_t1.txt")
rich_2016_t1 = readdlm("data/rich_mat_2016_t1.txt")
rich_2006_t2 = readdlm("data/rich_mat_2006_t2.txt")
rich_2016_t2 = readdlm("data/rich_mat_2016_t2.txt")

rich_2011_p_pi = readdlm("data/rich_mat_2011_p_pi.txt")
rich_2011_m_pi = readdlm("data/rich_mat_2011_m_pi.txt")
rich_2011_p_K = readdlm("data/rich_mat_2011_p_K.txt")
rich_2011_m_K = readdlm("data/rich_mat_2011_m_K.txt")
rich_2011_p_p = readdlm("data/rich_mat_2011_p_p.txt")
rich_2011_m_p = readdlm("data/rich_mat_2011_m_p.txt")

# 2006
pip2pip06 = [rich_2006_t1[:,3],rich_2006_t2[:,3]]
kp2pip06 = [rich_2006_t1[:,4],rich_2006_t2[:,4]]
pp2pip06 = [rich_2006_t1[:,5],rich_2006_t2[:,5]]
pip2kp06 = [rich_2006_t1[:,6],rich_2006_t2[:,6]]
kp2kp06 = [rich_2006_t1[:,7],rich_2006_t2[:,7]]
pp2kp06 = [rich_2006_t1[:,8],rich_2006_t2[:,8]]
pip2pp06 = [rich_2006_t1[:,9],rich_2006_t2[:,9]]
kp2pp06 = [rich_2006_t1[:,10],rich_2006_t2[:,10]]
pp2pp06 = [rich_2006_t1[:,11],rich_2006_t2[:,11]]
pim2pim06 = [rich_2006_t1[:,12],rich_2006_t2[:,12]]
km2pim06 = [rich_2006_t1[:,13],rich_2006_t2[:,13]]
pm2pim06 = [rich_2006_t1[:,14],rich_2006_t2[:,14]]
pim2km06 = [rich_2006_t1[:,15],rich_2006_t2[:,15]]
km2km06 = [rich_2006_t1[:,16],rich_2006_t2[:,16]]
pm2km06 = [rich_2006_t1[:,17],rich_2006_t2[:,17]]
pim2pm06 = [rich_2006_t1[:,18],rich_2006_t2[:,18]]
km2pm06 = [rich_2006_t1[:,19],rich_2006_t2[:,19]]
pm2pm06 = [rich_2006_t1[:,20],rich_2006_t2[:,20]]

# 2011
pip2pip11 = [rich_2011_p_pi[:,6],rich_2011_p_pi[:,10]]
kp2pip11 = [rich_2011_p_K[:,6],rich_2011_p_K[:,10]]
pp2pip11 = [rich_2011_p_p[:,6],rich_2011_p_p[:,10]]
pip2kp11 = [rich_2011_p_pi[:,7],rich_2011_p_pi[:,11]]
kp2kp11 = [rich_2011_p_K[:,7],rich_2011_p_K[:,11]]
pp2kp11 = [rich_2011_p_p[:,7],rich_2011_p_p[:,11]]
pip2pp11 = [rich_2011_p_pi[:,8],rich_2011_p_pi[:,12]]
kp2pp11 = [rich_2011_p_K[:,8],rich_2011_p_K[:,12]]
pp2pp11 = [rich_2011_p_p[:,8],rich_2011_p_p[:,12]]
pim2pim11 = [rich_2011_m_pi[:,6],rich_2011_m_pi[:,10]]
km2pim11 = [rich_2011_m_K[:,6],rich_2011_m_K[:,10]]
pm2pim11 = [rich_2011_m_p[:,6],rich_2011_m_p[:,10]]
pim2km11 = [rich_2011_m_pi[:,7],rich_2011_m_pi[:,11]]
km2km11 = [rich_2011_m_K[:,7],rich_2011_m_K[:,11]]
pm2km11 = [rich_2011_m_p[:,7],rich_2011_m_p[:,11]]
pim2pm11 = [rich_2011_m_pi[:,8],rich_2011_m_pi[:,12]]
km2pm11 = [rich_2011_m_K[:,8],rich_2011_m_K[:,12]]
pm2pm11 = [rich_2011_m_p[:,8],rich_2011_m_p[:,12]]

# 2016
pim2pim16 = [rich_2016_t1[:,3],rich_2016_t2[:,3]]
pim2km16 = [rich_2016_t1[:,4],rich_2016_t2[:,4]]
pim2pm16 = [rich_2016_t1[:,5],rich_2016_t2[:,5]]
pip2pip16 = [rich_2016_t1[:,6],rich_2016_t2[:,6]]
pip2kp16 = [rich_2016_t1[:,7],rich_2016_t2[:,7]]
pip2pp16 = [rich_2016_t1[:,8],rich_2016_t2[:,8]]
km2pim16 = [rich_2016_t1[:,9],rich_2016_t2[:,9]]
km2km16 = [rich_2016_t1[:,10],rich_2016_t2[:,10]]
km2pm16 = [rich_2016_t1[:,11],rich_2016_t2[:,11]]
kp2pip16 = [rich_2016_t1[:,12],rich_2016_t2[:,12]]
kp2kp16 = [rich_2016_t1[:,13],rich_2016_t2[:,13]]
kp2pp16 = [rich_2016_t1[:,14],rich_2016_t2[:,14]]
pm2pim16 = [rich_2016_t1[:,15],rich_2016_t2[:,15]]
pm2km16 = [rich_2016_t1[:,16],rich_2016_t2[:,16]]
pm2pm16 = [rich_2016_t1[:,17],rich_2016_t2[:,17]]
pp2pip16 = [rich_2016_t1[:,18],rich_2016_t2[:,18]]
pp2kp16 = [rich_2016_t1[:,19],rich_2016_t2[:,19]]
pp2pp16 = [rich_2016_t1[:,20],rich_2016_t2[:,20]]


for i in 1:2
    t = string(theta[i], " < ", theta[i+1])
    # PI+
    plot(mom06, pip2pip06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\pi^+ \rightarrow \pi^+$",
                              label = "2006",
                              title = t)
    plot!(mom11, pip2pip11[i], lw=3,
                               label = "2011")
    plot!(mom16, pip2pip16[i], lw=3,
                              label = "2016")
    if i==1
        savefig("plots/t1/pip2pip.png")
    else
        savefig("plots/t2/pip2pip.png")
    end

    plot(mom06, pip2kp06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\pi^+ \rightarrow K^+$",
                              label = "2006",
                              title = t)
    plot!(mom11, pip2kp11[i], lw=3,
                             label = "2011")
    plot!(mom16, pip2kp16[i], lw=3,
                              label = "2016")
    if i==1
        savefig("plots/t1/pip2kp.png")
    else
        savefig("plots/t2/pip2kp.png")
    end

    plot(mom06, pip2pp06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\pi^+ \rightarrow p$",
                              label = "2006",
                              title = t)
    plot!(mom11, pip2pp11[i], lw=3,
                             label = "2011")
    plot!(mom16, pip2pp16[i], lw=3,
                              label = "2016")
    if i==1
        savefig("plots/t1/pip2pp.png")
    else
        savefig("plots/t2/pip2pp.png")
    end

    # PI-
    plot(mom06, pim2pim06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\pi^- \rightarrow \pi^-$",
                              label = "2006",
                              title = t)
    plot!(mom11, pim2pim11[i], lw=3,
                             label = "2011")
    plot!(mom16, pim2pim16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/pim2pim.png")
    else
      savefig("plots/t2/pim2pim.png")
    end

    plot(mom06, pim2km06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\pi^- \rightarrow K^-$",
                              label = "2006",
                              title = t)
    plot!(mom11, pim2km11[i], lw=3,
                             label = "2011")
    plot!(mom16, pim2km16[i], lw=3,
                              label = "2016")
    if i==1
        savefig("plots/t1/pim2km.png")
    else
        savefig("plots/t2/pim2km.png")
    end

    plot(mom06, pim2pm06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\pi^- \rightarrow \bar p$",
                              label = "2006",
                              title = t)
    plot!(mom11, pim2pm11[i], lw=3,
                             label = "2011")
    plot!(mom16, pim2pm16[i], lw=3,
                              label = "2016")
    if i==1
        savefig("plots/t1/pim2pm.png")
    else
        savefig("plots/t2/pim2pm.png")
    end

    # K+
    plot(mom06, kp2pip06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$K^+ \rightarrow \pi^+$",
                              label = "2006",
                              title = t)
    plot!(mom11, kp2pip11[i], lw=3,
                             label = "2011")
    plot!(mom16, kp2pip16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/kp2pip.png")
    else
      savefig("plots/t2/kp2pip.png")
    end

    plot(mom06, kp2kp06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$K^+ \rightarrow K^+$",
                              label = "2006",
                              title = t)
    plot!(mom11, kp2kp11[i], lw=3,
                             label = "2011")
    plot!(mom16, kp2kp16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/kp2kp.png")
    else
      savefig("plots/t2/kp2kp.png")
    end

    plot(mom06, kp2pp06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$K^+ \rightarrow p$",
                              label = "2006",
                              title = t)
    plot!(mom11, kp2pp11[i], lw=3,
                             label = "2011")
    plot!(mom16, kp2pp16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/kp2pp.png")
    else
      savefig("plots/t2/kp2pp.png")
    end

    # K-
    plot(mom06, km2pim06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$K^- \rightarrow \pi^-$",
                              label = "2006",
                              title = t)
    plot!(mom11, km2pim11[i], lw=3,
                             label = "2011")
    plot!(mom16, km2pim16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/km2pim.png")
    else
      savefig("plots/t2/km2pim.png")
    end

    plot(mom06, km2km06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$K^- \rightarrow K^-$",
                              label = "2006",
                              title = t)
    plot!(mom11, km2km11[i], lw=3,
                             label = "2011")
    plot!(mom16, km2km16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/km2km.png")
    else
      savefig("plots/t2/km2km.png")
    end

    plot(mom06, km2pm06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$K^- \rightarrow \bar p$",
                              label = "2006",
                              title = t)
    plot!(mom11, km2pm11[i], lw=3,
                             label = "2011")
    plot!(mom16, km2pm16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/km2pm.png")
    else
      savefig("plots/t2/km2pm.png")
    end

    # P
    plot(mom06, pp2pip06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$p \rightarrow \pi^+$",
                              label = "2006",
                              title = t)
    plot!(mom11, pp2pip11[i], lw=3,
                             label = "2011")
    plot!(mom16, pp2pip16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/pp2pip.png")
    else
      savefig("plots/t2/pp2pip.png")
    end

    plot(mom06, pp2kp06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$p \rightarrow K^+$",
                              label = "2006",
                              title = t)
    plot!(mom11, pp2kp11[i], lw=3,
                             label = "2011")
    plot!(mom16, pp2kp16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/pp2kp.png")
    else
      savefig("plots/t2/pp2kp.png")
    end

    plot(mom06, pp2pp06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$p \rightarrow p$",
                              label = "2006",
                              title = t)
    plot!(mom11, pp2pp11[i], lw=3,
                             label = "2011")
    plot!(mom16, pp2pp16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/pp2pp.png")
    else
      savefig("plots/t2/pp2pp.png")
    end

    # PBAR
    plot(mom06, pm2pim06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\bar p \rightarrow \pi^-$",
                              label = "2006",
                              title = t)
    plot!(mom11, pm2pim11[i], lw=3,
                             label = "2011")
    plot!(mom16, pm2pim16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/pm2pim.png")
    else
      savefig("plots/t2/pm2pim.png")
    end

    plot(mom06, pm2km06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\bar p \rightarrow K^-$",
                              label = "2006",
                              title = t)
    plot!(mom11, pm2km11[i], lw=3,
                             label = "2011")
    plot!(mom16, pm2km16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/pm2km.png")
    else
      savefig("plots/t2/pm2km.png")
    end

    plot(mom06, pm2pm06[i], lw=3,
                              xlims = (0,50),
                              ylims = (0,1.1),
                              xlabel = L"$p_h$",
                              ylabel = L"$\bar p \rightarrow \bar p$",
                              label = "2006",
                              title = t)
    plot!(mom11, pm2pm11[i], lw=3,
                             label = "2011")
    plot!(mom16, pm2pm16[i], lw=3,
                              label = "2016")
    if i==1
      savefig("plots/t1/pm2pm.png")
    else
      savefig("plots/t2/pm2pm.png")
    end
end
