


from A_config.n01_O_config_params import *
from A_config.n02_O_config_analysis_functions import *



#### TRIG DICT
ntrail_dict_allpatient = {   'LH018' : {'RB' : 1, 'HV' : 1, 'AoRoCo': 1, 'AoR-Co': 3, 'A+RoCc': 2, 'A+R-Cc': 3, 'A+RoC-': 2, 'A+R-C-': 3},
                        }

ntrail_dict_allpatient_pilote = {   'NS217' : {'AoRoCo': 2, 'AoR-Co': 3, 'A+RoCc': 2, 'A+R-Cc': 3, 'A+RoC-': 1, 'A+R-C-': 2},
                        }



trig_dict_pilote = {'jules' :

                    #pilote 1
                    # {'AoRoCo' : [int(1.02e4), int(2.893e4)], 'AoR-Co01' : [int(2.893e4), int(4.736e4)], 'AoR-Co02' : [int(4.736e4), int(6.558e4)], 
                    # 'A+RoCc' : [int(1.0651e5), int(1.2376e5)], 'A+R-Cc01' : [int(1.346e5), int(1.4383e5)], 'A+R-Cc02' : [int(1.4565e5), int(1.6347e5)], 
                    # 'A+RoC-' : [int(1.8166e5), int(1.9985e5)], 'A+R-C-01' : [int(2.141e5), int(2.3107e5)]}
        
                    #pilote 2
                    # {'AoRoCo' : [int(3.1e2), int(1.306e4)], 'AoR-Co01' : [int(1.306e4), int(2.592e4)], 'AoR-Co02' : [int(2.675e4), int(3.904e4)], 
                    # 'A+RoCc' : [int(4.156e4), int(5.361e4)], 'A+R-Cc01' : [int(5.638e4), int(6.837e4)], 'A+R-Cc02' : [int(6.837e4), int(8.131e4)], 
                    # 'A+RoC-' : [int(1.1459e5), int(1.2626e5)], 'A+R-C-01' : [int(1.2801e5), int(1.4e5)], 'A+R-C-02' : [int(1.4064e5), int(1.5306e5)]}
        
                    #pilote 3
                    {'AoRoCo' : [int(8404), int(26840)], 'AoR-Co01' : [int(28186), int(46932)], 'AoR-Co02' : [int(52238), int(70849)], 
                    'A+RoCc' : [int(147106), int(165534)], 'A+R-Cc01' : [int(167138), int(185370)], 'A+R-Cc02' : [int(186461), int(204680)], 
                    'A+RoC-' : [int(207969), int(226702)], 'A+R-C-01' : [int(227400), int(245644)], 'A+R-C-02' : [int(246424), int(265340)]},
        
                    'jose' :
        
                    #pilote 1
                    {'AoRoCo' : [int(41403), int(59933)], 'AoR-Co01' : [int(114820), int(133687)], 'AoR-Co02' : [int(135718), int(154151)], 
                    'A+RoCc' : [int(155936), int(174080)], 'A+R-Cc01' : [int(175920), int(194345)], 'A+R-Cc02' : [int(195225), int(215060)], 
                    'A+RoC-' : [int(267012), int(285473)], 'A+R-C-01' : [int(287389), int(306136)], 'A+R-C-02' : [int(307220), int(325884)]},

                    'test_aaron_params' :
        
                    #pilote 1
                    {'AoRoCo' : [int(0), int(59000)]},

                    'nose_mouth' :

                    {'new_nose' : {'AoRo' : [int(9.5e2), int(1.934e4)], 'AoR-' : [int(2.076e4), int(3.904e4)], 'A+Ro' : [int(4.011e4), int(5.834e4)]}, 
        
                    'old_nose' : {'AoRo' : [int(6.275e4), int(8.086e4)], 'AoR-' : [int(8.172e4), int(1.007e5)], 'A+Ro' : [int(1.0132e5), int(1.201e5)]},
        
                    'mouth' : {'AoRo' : [int(1.3516e5), int(1.5410e5)], 'AoR-' : [int(1.5581e5), int(1.7426e5)], 'A+Ro' : [int(1.7602e5), int(1.946e5)]}},
        
                    }

################################
######## AUXCHAN ########
################################

auxchan_allsujet = {

    'LH018' : ['TRIG', 'OSAT', 'PR', 'Pleth', 'respi', 'O2', 'CO2'],



}







