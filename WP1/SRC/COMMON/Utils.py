TESTING = True
'''
Be smart, don't be a retard!

# Ver satélite 30
gawk '$5 == 30' OBS_TLSZ_Y19D014.dat

# Ver satélite 30 e imprimir líneas 1 a 10
gawk  '$5 == 30' OBS_TLSZ_Y19D014.dat | sed -n '1,10p;11q'
'''
def get_sat_label(PRN):
    return "G" + "%02d" % int(PRN)

def set_sat_valid(a_sat='', a_activate=False, a_rejection_cause=0, PreproObsInfo=dict()):
    if a_activate:
        # Method only for deactivate a satellite =)
        return

    PreproObsInfo[a_sat]["ValidL1"] = 0
    PreproObsInfo[a_sat]["RejectionCause"] = a_rejection_cause
