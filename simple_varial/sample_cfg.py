# Configuration file for 2018 samples

# The name of the TTree in the ntuple.
treename = 'ntuple'

# The stacking order of physics processes in the plots.
# The list is ordered from top to bottom of the stack.
stacking_order = [
'FakeMC',
'FakeSig',
]

# The fill and line colors assigned to each physics process
# in the plots. The current palette copies HIG-16-044.
sample_colors = {
'FakeMC' : 596,
'FakeSig' : 632,
}

# The input pattern used to glob for sample files given their input token.
input_pattern = ['/data/pubfs/zhanglic/workspace/HTauTauBSM/2018/*/mt/*%s*.root', '/data/pubfs/zhanglic/workspace/HTauTauBSM/2018/*/mt/*%s*.root']



def get_samples(channel, signal_overlay=True, **kwargs):
    
    sf_lumi = 1.
    # Samples common to all channels.
    samples = {
        'Data_18_test_MC': ['4==4', 1., 'FakeMC', ['ElectronEmbedding_Run2018A','ElectronEmbedding_Run2018B','ElectronEmbedding_Run2018C','ElectronEmbedding_Run2018D'],0],
    }

    # Samples specific to the Zll channel.

    if channel == 'Zll':
        samples.update({
            'Data_18_test': ['4==4', 1., 'Data', ['ElectronEmbedding_Run2018A','ElectronEmbedding_Run2018B','ElectronEmbedding_Run2018C','ElectronEmbedding_Run2018D'],1],
        })

        if signal_overlay:
            samples.update({
                'Data_18_test_sig': ['4==4', 0.5, 'FakeSig', ['ElectronEmbedding_Run2018A','ElectronEmbedding_Run2018B','ElectronEmbedding_Run2018C','ElectronEmbedding_Run2018D'],1],
            })

    return samples

