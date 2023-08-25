# functions of HOTRG

using LinearAlgebra
function getMPO(g,bt,N)        
    sx = [0 1; 1 0]; sy = [0 -im; im 0]; sz = [1 0; 0 -1]; sI = [1 0; 0 1];
    hlocal = -(kron(sz,sz) + g*kron(sx,sI));

    # exponentiate the hamilton
    d = 2;
    gateAB = reshape(exp(-bt/N*hlocal),d,d,d,d);
    gateBA = reshape(exp(-bt/N*hlocal),d,d,d,d);

    # transform gate to local MPO: TA, TB
    gateABp = permutedims(gateAB,[1,3,2,4]);
    AB = svd(reshape(gateABp,d^2,d^2));
    gateBAp = permutedims(gateBA,[1,3,2,4]);
    BA = svd(reshape(gateBAp,d^2,d^2));

    ABa = reshape(AB.U * sqrt.(diagm(0 => AB.S)),  d,d,d^2);
    ABb = reshape(sqrt.(diagm(0 => AB.S)) * AB.Vt, d^2,d,d);
    BAa = reshape(sqrt.(diagm(0 => AB.S)) * BA.Vt, d^2,d,d);
    BAb = reshape(BA.U * sqrt.(diagm(0 => BA.S)),  d,d,d^2);

    TA = ncon(Any[BAa,ABa], Any[[-1,-3,1],[1,-4,-2]]);
    TB = ncon(Any[ABb,BAb], Any[[-1,1,-4],[-3,1,-2]]);

    # contract TA,TB to a two-sites MPO
    Tp = ncon(Any[TA,TB], Any[[-1,1,-3,-5],[1,-2,-4,-6]]);
    dT = size(Tp);
    T = reshape(Tp, dT[1],dT[2],dT[3]*dT[4],dT[5]*dT[6]);
    return T
end

function getGauge(T,chi)
    # do HOSVD to look for UL
    dT = size(T);
    MM = ncon(Any[T,T], Any[[-1,-3,-5,1],[-2,-4,1,-6]]);
    dM = size(MM);
    ML = reshape(MM,dM[1]*dM[2],dM[3]*dM[4]*dM[5]*dM[6]);
    ML = ML * ML';
    ML = 0.5* real(ML + ML');
    L = svd(ML);
    errL = norm(L.S[chi+1:end]);
    SL = diagm(0 => L.S);
    
    # do HOSVD to look for UR
    MRp = permutedims(MM,[3,4,1,2,5,6]);
    MR = reshape(MRp,dM[3]*dM[4],dM[1]*dM[2]*dM[5]*dM[6]);
    MR = MR*MR';
    MR = 0.5* real(MR + MR');
    R = svd(MR);
    errR = norm(R.S[chi+1:end]);
    SR = diagm(0 => R.S);
    
    # compare the error of truncation between UL and UR to decide gauge 
    if errL < errR
        err = errL;
        gauge = L.U[:, setdiff(1:end, chi+1:end)]; # delete multiple columns
        chi = min(chi,dM[1]*dM[2]);
        gauge = reshape(gauge, dM[1],dM[2],chi);
    else
        err = errR;
        gauge = R.U[:, setdiff(1:end, chi+1:end)]; 
        chi = min(chi,dM[3]*dM[4]);
        gauge = reshape(gauge,dM[3],dM[4],chi);
    end
    return gauge
end

# calculate partition function of a 1D-transverse-Ising-Model
function coarseGrain(g,bt,N,chi,L)
    # g: transverse field
    # bt: beta
    # N: steps of TEBD(convert to N rows of HOTRG)
    # chi: low rank of gauge
    # L: L cells of local AB 
        n = Int(log2(N));
        T = getMPO(g,bt,N);
        for i = 1:n
            gauge = getGauge(T,chi);
            T = ncon(Any[gauge,T,T,gauge], 
                    Any[[1,3,-1],[1,2,-3,5],[3,4,5,-4],[2,4,-2]]);
        end
        # partial trace of last tensor
        p = ncon(Any[T],Any[-1,-2,1,1]);
        Z = tr(p^L);
        return Z
    end