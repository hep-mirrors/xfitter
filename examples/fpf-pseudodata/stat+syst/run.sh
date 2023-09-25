# The commands to run Hessian PDF profiling using FPF pseudodata
# Using estimates for statistical and systematic uncertainties here.
# Running the stat-unc-only case in the "stat" directory is prerequisite for
# producing the comparison plot in this script.

ln -s ../../../datafiles

DETECTOR="FPF"
DETSIMPLE="FPF"
UNCLVL="stat+syst"

PDF="PDF4LHC21"  #N.B. filenames only! Change also in parameters.yaml
PDFLONG="PDF4LHC21_40"

PLOTOPTS="--splitplots-pdf --therr --bands --q2all --relative-errors --xrange 1e-3:0.6"

# Run xFitter:
xfitter

# Visualize results:
mv plots/plots.pdf plots/plots_OLD.pdf
xfitter-draw ${PLOTOPTS} output:"Baseline (BL)" profile:../stat/output:"BL+${DETECTOR}, stat" profile:output:"BL+${DETECTOR}, stat+syst"
mv plots/plots.pdf plots/${PDF}_vs_profiled.pdf

cp plots/q2_10000_pdf_uv_ratio.pdf  ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_uv_ratio.pdf 
cp plots/q2_10000_pdf_dv_ratio.pdf  ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_dv_ratio.pdf 
cp plots/q2_10000_pdf_g_ratio.pdf   ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_g_ratio.pdf
cp plots/q2_10000_pdf_Sea_ratio.pdf ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_Sea_ratio.pdf
cp plots/q2_10000_pdf_s_ratio.pdf   ${UNCLVL}_${DETSIMPLE}_q2_10000_pdf_s_ratio.pdf

# Produce a new PDF in LHAPDF6 format:
#xfitter-process profile output/pdf_shifts.dat output/pdf_rotate.dat `lhapdf-config --datadir`/${PDFLONG}/ ${PDF}-${DETECTOR}

# Save the new PDF set into LHAPDF6 collection:
#cp -r ${PDF}-profiled/ `lhapdf-config --datadir`
