#!/bin/csh -f

set CUT = -6
set tempselect = 3500.
set weakratio = 1.e${CUT}

set LAMBDAMIN = '7000'
set LAMBDAMAX = '9000'

set DIRECTORY = '7000-9000'

mkdir ${DIRECTORY}


#------------------------------------------------------------------
# Selection de raies CN pour Bsyn. BPz 18/02-99
#------------------------------------------------------------------

translatelinelists_autocount_identif <<eof
linelist/linelistCN_alliso_R240511.dat
C2N4sg
${DIRECTORY}/C12N14R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
linelist/linelistCN_alliso_R240511.dat
C2N5sg
${DIRECTORY}/C12N15R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
linelist/linelistCN_alliso_R240511.dat
C3N4sg
${DIRECTORY}/C13N14R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
linelist/linelistCN_alliso_R240511.dat
C3N5sg
${DIRECTORY}/C13N15R_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
linelist/linelistCN1214V130710.dat
C12N14
${DIRECTORY}/C12N14V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
linelist/linelistCN1215V130710.dat
C12N15
${DIRECTORY}/C12N15V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
linelist/linelistCN1314V130710.dat
C13N14
${DIRECTORY}/C13N14V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof

translatelinelists_autocount_identif <<eof
linelist/linelistCN1315V130710.dat
C13N15
${DIRECTORY}/C13N15V130710_${LAMBDAMIN}-${LAMBDAMAX}.list
${LAMBDAMIN}
${LAMBDAMAX}
${tempselect}
${weakratio}
eof


