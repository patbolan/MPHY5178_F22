%function DM = blochRK4_loop( N, B1TX, B1TY, Tp, OFFSET, T1, T2, DMO)
function Mfinal = blochRK4_loop(B1TX, B1TY, Tp, OFFSET, R1, R2, Minit)


Npts = max(size(B1TX(:)));
deltaT = Tp / (Npts+1);

Mnext = zeros(Npts, 3);
Mnext(1,:) = Minit;
for ndx=2:Npts
    Mnext(ndx,:) = blochRK4(Mnext(ndx-1,:), B1TX(ndx), B1TY(ndx), OFFSET, R1, R2, deltaT);
end

Mfinal = squeeze(Mnext(Npts,:));