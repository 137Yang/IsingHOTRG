using Plots, LaTeXStrings
@load "F:\\julia\\IsingHOTRG\\figure\\Z.jld2" btR btI Z
    # gr();
    plot_size = (3500,2625);
    contour(btR, btI, real(Z), levels=[0], seriescolor=:reds, legend=:none, size=plot_size)
    contour!(btR, btI, imag(Z), levels=[0], seriescolor=:blues, legend=:none, size=plot_size)
    xlabel!(L"β_r ")
    ylabel!(L"β_i(π)")
    #save plot as PNG
    savefig("F:\\julia\\IsingHOTRG\\figure\\z.png")