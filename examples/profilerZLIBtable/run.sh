# Perform a PDF Hessian profiling

PDF="CT14nnlo"
#PLOTOPTS="--asym --therr --bands --q2all --relative-errors"
PLOTOPTS="--splitplots-eps --therr --bands --q2all --relative-errors --scale68 --xrange 1e-4:0.7"

# Run xFitter:
xfitter

# Visualize results:
mv plots/plots.pdf plots/plots_OLD.pdf
xfitter-draw ${PLOTOPTS} output:${PDF} profile:output:jets-${PDF}-profiled
mv plots/plots.pdf plots/${PDF}_NNLO_vs_jets-profiled.pdf

# Produce a new PDF in lhapdf6 format:
#xfitter-process profile output/pdf_shifts.dat output/pdf_rotate.dat `lhapdf-config --datadir`/${PDF}/ ${PDF}-profiled

# Save the new PDF set into our lhapdf6 collection:
#cp -r ${PDF}-profiled/ `lhapdf-config --datadir`

