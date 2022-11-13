from sample_cfg import *
import sys
import varial
import numpy as np

####################
# General Settings #
####################

varial.settings.max_num_processes = 5
varial.settings.rootfile_postfixes += ['.pdf', '.png']
varial.settings.stacking_order = stacking_order
# try to find output on disk and don't run a step if present
enable_reuse_step = True

#################
# Plot Settings #
#################

name = 'test'

# weight = 'weight*DL_PreFitSF_4b'
weight = 'genWeight'

plot_vars = {    
    'mt_1'       :    ('mt_1',     ';m_{T}^{#mu} [GeV];', 28,0,140),
}

# from jet
plot_vars.update({
})

#######################################
# Samples, Selections, and Categories #
#######################################
btag = "( nbtag >= 1 ) "
tight_mT = " ( mt_1 < 40.0) "
opposite_sign = ' ((q_1 * q_2) < 0) '
tau_selections = " (id_tau_vsMu_Tight_2 > 0 && id_tau_vsJet_Medium_2 > 0 &&  id_tau_vsEle_VVLoose_2 > 0 ) "

the_samples_dict = get_samples(
    channel='Zll',
    signal_overlay=True,
)
regions = {
    "btag_tight_mT"         : '{0} && {1} && {2} && {3}'.format(btag, tight_mT, opposite_sign, tau_selections),
}

lepton_veto = "extramuon_veto == 0  && extraelec_veto == 0 " 
muon_selections = " pt_1> 25.0 "

selections = [
    lepton_veto, muon_selections # DR selections
]

the_category_dict = {
    'Htautau': [regions, selections, plot_vars],
}

# Blinding the Signal Region
def additional_input_hook(wrps):

    @varial.history.track_history
    def rebin_custom(w):        
        if w.in_file_path.startswith('btag_tight_mT'):
            if 'tau_pt' in w.name:
                print('REBIN Datacards in %s' % w.in_file_path)
                w.histo = w.histo.Rebin(2,'new name',[0,1,2])
        return w
    
    wrps = (rebin_custom(w) for w in wrps)
    return wrps
