TESTING = True
'''
Be smart, don't be a retard!

# Ver satélite 30
gawk '$5 == 30' OBS_TLSZ_Y19D014.dat

# Ver satélite 30 e imprimir líneas 1 a 10
gawk  '$5 == 30' OBS_TLSZ_Y19D014.dat | sed -n '1,10p;11q'
gawk  '$5 == 30 && $1>1000 && $1<3600 ' OBS_TLSZ_Y19D014.dat 

# Sacar la columna SOD y C1SMOOTHED para PRN10 y compararlas con la solucion
gawk '$4==10' PREPRO_OBS_TLSZ_Y19D014.dat | cut -c1-6 > hgm_sod_prn10.dat
gawk '$4==10' PREPRO_OBS_TLSZ_Y19D014.dat | cut -c66-82 > hgm_smoothedc1_prn10.dat
paste hgm_sod_prn10.dat hgm_smoothedc1_prn10.dat > hgm_prn10_plot_c1smoothed.dat
# SOLUTION
gawk '$4==10' SOLUTION.dat | cut -c1-6 > solution_sod_prn10.dat
gawk '$4==10' SOLUTION.dat | cut -c66-82 > solution_smoothedc1_prn10.dat
paste solution_sod_prn10.dat solution_smoothedc1_prn10.dat > solution_prn10_plot_c1smoothed.dat
# JOIN FILES
paste hgm_prn10_plot_c1smoothed.dat solution_smoothedc1_prn10.dat > C1SMOOTHED.dat
# PLOT
cat << EOF | gnuplot --persist
pipe heredoc> plot "C1SMOOTHED.dat" u 1:2, "" u 1:3
pipe heredoc> EOF

# Sacar columna SOD
cut -c1-6 INP.dat > OUT.dat
# Sacar columna C1SMOOTHED
cut -c66-82 SOLUTION.dat > solution_c1smoothed.dat

# Ver gaps
gawk '$8==6' SOLUTION.dat 
# Ver CS
gawk '$8==5' SOLUTION.dat 

# Para poder correr un script recién creado
chmod u+x plot_c1smoothed.sh 
'''
def get_sat_label(PRN):
    return "G" + "%02d" % int(PRN)

def set_sat_valid(PreproObs, a_valid, a_rejection_cause, ):
    PreproObs["ValidL1"] = a_valid
    PreproObs["RejectionCause"] = a_rejection_cause
