# 110917
# script obtained from https://www.biostars.org/p/114183/
# script is for filtering out by TLEN = signed observed Template LENgth. If all segments are mapped to the same reference, the unsigned observed template length equals the number of bases from the leftmost mapped base to the rightmost mapped base. 
# effectively, this filters paired end reads by fragment length, I set this as between 120bp and 170bp i.e. 150bp +/- 20% (mononucleosome = ~147bp)
# SIZE1 is the minimum template lenght
BEGIN { FS="\t"; SIZE1=170; S1=SIZE1*SIZE1; SIZE2=120; S2=SIZE2*SIZE2 }

# Write the headers 
/^@/ { print $0; next }

# Print only mononucleosomal entries
{ if ($9*$9 < S1 && $9*$9 > S2) print $0}
