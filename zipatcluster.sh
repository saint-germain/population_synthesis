NAME=$1
cp $NAME.e* results/
cp $NAME.o* results/
zip $NAME.zip results/res* results/*.e* results/*.o* results/info* results/header.txt parameters.in mypy.pbs main_ccx_par.py input.py
