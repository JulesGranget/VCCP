


from n00_config_params import *
from n00bis_config_analysis_functions import *


#### BLOCK ORDER
patient_block_order = {'NS217' : ['AoRoCo', 'A+RoCo', 
                                  'AoRoCo', 'AoR-Co', 'AoR-Co', 'AoR-Co',
                                  'A+RoCc', 'A+R-Cc', 'A+R-Cc', 'A+RoCc', 'A+R-Cc',
                                  'A+RoC-', 'A+R-C-', 'A+R-C-']}








#### TRIG DICT
trig_dict_allpatient = {   'NS217' : {
                'AoRoCo01' : [int(106212), int(118539)], 'A+RoCo00' : [int(130276), int(140631)], 
                'AoRoCo02' : [int(148047), int(160182)], 'AoR-Co01' : [int(215542), int(223640)], 'AoR-Co02' : [int(234101), int(246452)], 'AoR-Co03' : [int(252129), int(264546)], 
                'A+RoCc01' : [int(270372), int(282718)], 'A+R-Cc01' : [int(291832), int(304402)], 'A+R-Cc02' : [int(336897), int(348886)], 'A+RoCc02' : [int(357287), int(369303)], 'A+R-Cc03' : [int(382814), int(394845)],
                'A+RoC-01' : [int(415978), int(428044)], 'A+R-C-01' : [int(436847), int(448950)], 'A+R-C-02' : [int(460005), int(472243)]}
                
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
                    'A+RoC-' : [int(267012), int(285473)], 'A+R-C-01' : [int(287389), int(306136)], 'A+R-C-02' : [int(307220), int(325884)]}
        
                    }








#### DATA FILE NAME
datafile_name_allsujet = {'NS217' : {'B' : 'FirstDeriv_Neural-251110-120653_VCCP_NS217_B30_B3', 'P' : 'FirstDeriv_Neural-251110-120653_VCCP_NS217_B30_P1'}}








