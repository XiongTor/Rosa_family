#!/usr/bin/env Rscript
# Author: #!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-10-24
# Description: To simlulate different

# ==== Main Script Start ====
library(phybase)

# mptree1 = "(Gbi:8.4,(Atr:7.4,(Efe:7.2,((Peq:4.0,(Mac:3.4,Osa:3.4):0.55):2.1,((Lch:4.6,(Pam:1,Cka:1):3.6):1.3,(Cde:5.9,((Vvi:3.1,(Ath:2.7,Ppe:2.7):0.40):1.7,Aco:4.8):1.0):0.09):0.08):1.07):0.21):1);"
mptree1 = "(((((((((((((Sibbaldia_parviflora:0.012744,Sibbaldianthe_bifurca:0.019597)1.000000:0.011331,Alchemilla_faeroensis:0.076049)0.994427:0.001734,Comarum_palustre:0.015885)0.934500:0.010986,Fragaria_nilgerrensis:0.016319)0.987878:0.001493,(((Drymocallis_arguta:0.007451,Chamaecallis_perpusilloides:0.005753)1.000000:0.006421,Chamaerhodos_erecta:0.033023)0.950221:0.001202,(Dasiphora_fruticosa:0.005309,Potaninia_mongolica:0.00884)1.000000:0.012941)1.000000:0.013618)1.000000:0.0177,(Argentina_anserina:0.039105,Potentilla_acaulis:0.028777)1.000000:0.01601)1.000000:0.013833,((((Polylepis_tarapacana:0.006722,Acaena_ovalifolia:0.004862)0.997993:0.003289,Margyricarpus_pinnatus:0.005207)1.000000:0.013052,(Sanguisorba_minor:0.002133,Bencomia_exstipulata:0.000314)1.000000:0.024044)1.000000:0.02154,Agrimonia_pilosa_var._pilosa:0.044563)1.000000:0.018667)0.948422:0.00137,Rosa_chinensis:0.01352)1.000000:0.022532,((Coluria_longifolia:0.021155,Geum_urbanum:0.002677)0.585990:0.001711,(Taihangia_rupestris:0.001545,Waldsteinia_ternata:0.006699)1.000000:0.023664)1.000000:0.069627)0.999989:0.016996,Rubus_argutus:0.028827)1.000000:0.014548,Filipendula_ulmaria:0.073819)1.000000:0.041494,((((Cercocarpus_ledifolius:0.005756,Purshia_tridentata:0.009686)0.970670:0.006862,Chamaebatia_foliolosa:0.013455)0.999158:0.00432,Dryas_ajanensis_subsp._ajanensis:0.001308)1.000000:0.036657,((Neillia_sinensis:0.010776,Physocarpus_opulifolius:0.003152)1.000000:0.035451,((Lyonothamnus_floribundus:0.05465,Prunus_padus:0.034773)0.967235:0.003392,((Aruncus_dioicus:0.019662,Spiraea_blumei:0.030565)1.000000:0.044091,(((((Prinsepia_uniflora:0.01307,Exochorda_racemosa_subsp._serratifolia:0.000128)0.694596:0.020665,Oemleria_cerasiformis:0.002035)1.000000:0.04866,((Coleogyne_ramosissima:0.016216,Kerria_japonica:0.015361)1.000000:0.015848,Rhodotypos_scandens:0.017639)1.000000:0.01436)0.999655:0.003648,Sorbaria_kirilowii_var._arborea:0.061299)0.995053:0.00242,(Gillenia_trifoliata:0.017443,(((Lindleya_mespiloides:0.011519,Vauquelinia_australis:0.005602)0.406233:0.022501,Kageneckia_oblonga:0.027559)0.401694:0.000534,(((Hesperomeles_cuneata:0.013523,Phippsiomeles:0.02771)0.673232:0.002402,Dichotomanthes_tristaniicarpa:0.00589)0.403753:0.019933,(Cormus_domestica:0.027516,(Chamaemeles_coriacea:0.023068,(Cotoneaster_frigidus:0.018199,(Stranvaesia_nussia:0.007368,((((Eriobotrya_japonica:0.00013,Rhaphiolepis_ferruginea:0.005151)0.999900:0.006838,(Sorbus_aucuparia:0.010966,((Pyrus_communis:0.014678,(((Karpatiosorbus_bristoliensis:0.014797,Aria_edulis:0.005206)0.742507:0.000955,Torminalis_glaberrima:0.005611)0.627254:0.0044,((Thomsonaria_caloneura:0.011519,Griffitharia_hemsleyi:0.012004)0.379936:0.000568,(Micromeles_alnifolia:0.003272,Alniaria_alnifolia:0.003898)0.999986:0.005591)0.795737:0.005853)0.982549:0.002308)0.617063:0.000726,(((((Aronia_arbutifolia:0.007572,Osteomeles_schweriniae:0.017053)0.521730:0.000518,(Pseudocydonia_sinensis:0.008932,Cydonia_oblonga:0.006363)0.999992:0.003439)0.607332:0.008439,Chaenomeles_speciosa:0.015018)0.422989:0.00891,Pourthiaea_amphidoxa:0.021504)0.539493:0.000688,Pyracantha_coccinea:0.011258)0.410354:0.001751)0.469162:0.004069)0.554703:0.002064)0.416980:0.000538,Photinia_prunifolia:0.004825)0.546728:0.017268,((((Malacomeles_denticulata:0.005189,Amelanchier_laevis:0.006411)0.999958:0.012534,Peraphyllum_ramosissimum:0.005455)0.806890:0.001022,Crataegus_hupehensis:0.01584)0.995180:0.002215,((Macromeles_tschonoskii:0.00438,Malus_domestica:0.004435)1.000000:0.001049,Weniomeles_bodinieri:0.004552)0.924394:0.016055)0.359473:0.005166)0.467068:0.000647)0.515927:0.010251)0.546501:0.00091)0.514927:0.000937)0.542910:0.001218)0.999995:0.009996)1.000000:0.009443)1.000000:0.044049)0.723670:0.004036)0.868523:0.001268)0.996614:0.005311)1.000000:0.023904)0.569203:0.001571)1.000000:0.038285,(Elaeagnus_angustifolia:0.175192,(Zelkova_schneideriana:0.075722,Morus_indica:0.088154)1.000000:0.02664)1.000000:0.038285);"


spname <- species.name(mptree1)
nodematrix <- read.tree.nodes(str=mptree1, name=spname)$nodes

theta_levels <- c(0.001,0.055,0.840,8.845)   # 不同程度的 ILS（从弱到强）
n_gene <- 2000                    # 每组生成2000棵基因树

for (theta in theta_levels) {
  nodematrix[,5] <- theta
  
  genetrees <- character(n_gene)
  for (i in 1:n_gene) {
    genetrees[i] <- sim.coaltree.sp(rootnode=max(nodematrix[,1]),
                                    nodematrix=nodematrix,
                                    nspecies=length(spname),
                                    seq=rep(1, length(spname)),
                                    name=spname)$gt
  }
  
  outname <- paste0("sim_theta_", theta, ".tre")
  write(genetrees, outname)
  cat("Done: theta =", theta, "\n")
}

# Date: 2025-10-30
# Description: 

# ==== Main Script Start ====

