function logT = f_logT(e_VG,iGraf)

persistent logTloc

if e_VG.istep==0
    logT = 0;
    logTloc = 0;
else
    if iGraf==1
        logT = logTloc + e_VG.Dtime;
        logTloc = logT;
    else
        logT = logTloc;
    end
end