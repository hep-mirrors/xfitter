# Contains all the commands needed for profiling the Contact Interactions (CI)
# Wilson coefficient c_1, using the left-handed CI model (CIcase = 1) and new
# physics scale 10 TeV (LamoTeV = 10)

NAME="CT14nlo_CI1L10"

#Need --scale68 for prodiling with CT PDFs
PLOTOPTS="--therr --bands --q2all --relative-errors --scale68"

# Run xFitter:
xfitter

# Visualize results:
mv plots/plots.pdf plots/plots_OLD.pdf
xfitter-draw ${PLOTOPTS} output:${NAME} profile:output:${NAME}-profiled

# Produce a new PDF in lhapdf6 format:
#xfitter-process profile output/pdf_shifts.dat output/pdf_rotate.dat `lhapdf-config --datadir`/${NAME}/ ${NAME}-profiled

# Save the new PDF set into our lhapdf6 collection:
#cp -r ${NAME}-profiled/ `lhapdf-config --datadir`
