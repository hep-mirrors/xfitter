cat <<END > head.LHgrid
'Version' '5.8'
'Description: '
'HERAPDF10 fit: varaint'
'member 0 - central'
' '
'Alphas:'
'Variable','nlo','EvolCode'
1,91.187,1.40,4.75,180.0
'MinMax:'
20,1
1.E-06,1.,1.0,200000000.
'QCDparams:'
0,1
0.338,0.243
'Parameterlist:'
'list',0,1
0.1176
'Evolution:'
'nlo',1.9,1.0
'HERAGRID10'
0,0
END

echo "Parse..."
bin/tolhapdf
cat head.LHgrid lhapdf_tail.dat > PDFs.LHgrid
echo "'End:'" >> PDFs.LHgrid

rm -f head.LHgrid lhapdf_tail.dat

ls -l PDFs.LHgrid


