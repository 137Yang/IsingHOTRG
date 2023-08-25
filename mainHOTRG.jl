include("parinput.jl")
include("doHOTRG.jl")
include("ncon.jl")
using JLD2, Plots, LaTeXStrings
    gi, btR, btI, N, chi, L = parinput()
    Ng = length(gi);
    Nbr = length(btR);
    Nbi = length(btI)
    
        Z = zeros(ComplexF64, Nbi,Nbr);
        g = gi
        for i in 1:Nbi
            for j in 1:Nbr
                x = btR[j];
                y = btI[i];
                bt = x + y*π*im;
                Z[i,j] = coarseGrain(g,bt,N,chi,L);
                @save "F:\\julia\\IsingHOTRG\\figure\\Z.jld2" btR btI Z
            end
        end       
    @load "F:\\julia\\IsingHOTRG\\figure\\Z.jld2" btR btI Z
    # gr();
    plot_size = (3500,2625);
    contour(btR, btI, real(Z), levels=[0], seriescolor=:reds, legend=:none, size=plot_size)
    contour!(btR, btI, imag(Z), levels=[0], seriescolor=:blues, legend=:none, size=plot_size)
    xlabel!(L"β_r ")
    ylabel!(L"β_i(π)")
    #save plot as PNG
    savefig("F:\\julia\\IsingHOTRG\\figure\\z.png")