function parinput()
    btR = range(0, 1, length=50);
    btI = range(-1, 1, length=100);
    chi = 5;
    nstep = 2;
    N = 2^nstep;
    nsite = 4;
    L = Int(nsite/2);
    gi = 0.1;
    return gi, btR, btI, N, chi, L 
end
