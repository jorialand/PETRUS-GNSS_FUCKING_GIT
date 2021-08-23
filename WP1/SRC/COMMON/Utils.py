def get_sat_label(PRN):
    return "G" + "%02d" % int(PRN)

def set_sat_valid(a_sat='', a_activate=False, a_rejection_cause=0, PreproObsInfo=dict()):
    if a_activate:
        # Method only for deactivate a satellite =)
        return

    PreproObsInfo[a_sat]["ValidL1"] = 0
    PreproObsInfo[a_sat]["RejectionCause"] = a_rejection_cause
