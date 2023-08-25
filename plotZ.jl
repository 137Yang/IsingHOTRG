using Plots, LaTeXStrings
contour(btR, btI, real(Z), levels=[0], color=[:red], legend=:none)
contour!(btR, btI, imag(Z), levels=[0], color=[:blue], legend=:none)
xlabel!(L"β_r ")
ylabel!(L"β_i(π)")