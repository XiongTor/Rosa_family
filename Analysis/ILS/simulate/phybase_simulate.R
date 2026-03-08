#!/usr/bin/env Rscript
# Author: #!/usr/bin/env Rscript
# Author: Tao Xiong
# Date: 2025-10-24
# Description: To simlulate different

# ==== Main Script Start ====
library(phybase)

# mptree1 = "(Gbi:8.4,(Atr:7.4,(Efe:7.2,((Peq:4.0,(Mac:3.4,Osa:3.4):0.55):2.1,((Lch:4.6,(Pam:1,Cka:1):3.6):1.3,(Cde:5.9,((Vvi:3.1,(Ath:2.7,Ppe:2.7):0.40):1.7,Aco:4.8):1.0):0.09):0.08):1.07):0.21):1);"
mptree1 = "((((((((((((Sibbaldia_parviflora:0.021117,Sibbaldianthe_bifurca:0.032293)1.000000:0.00446,Alchemilla_faeroensis:0.086909)1.000000:0.002667,Comarum_palustre:0.02496)0.999966:0.003304,Fragaria_nilgerrensis:0.01589)0.998665:0.001845,(((Chamaecallis_perpusilloides:0.007008,Drymocallis_arguta:0.009677)1.000000:0.004292,Chamaerhodos_erecta:0.038027)0.999995:0.000415,(Dasiphora_fruticosa:0.007339,Potaninia_mongolica:0.015046)1.000000:0.007818)1.000000:0.007427)1.000000:0.007169,(Potentilla_acaulis:0.037517,Argentina_anserina:0.041674)1.000000:0.005684)1.000000:0.000468,(((((Margyricarpus_pinnatus:0.009446,Acaena_ovalifolia:0.002285)0.925528:0.001274,Polylepis_tarapacana:0.004565)1.000000:0.010586,(Bencomia_exstipulata:0.004653,Sanguisorba_minor:0.001838)1.000000:0.01295)1.000000:0.011806,Agrimonia_pilosa_var._pilosa:0.04894)1.000000:0.006456,Rosa_chinensis:0.019311)0.844168:0.001716)0.963873:0.014239,((Coluria_longifolia:0.018569,Geum_urbanum:0.030186)0.527904:0.004713,(Waldsteinia_ternata:0.010945,Taihangia_rupestris:0.0035)0.645983:0.005365)0.714809:0.02059)0.612414:0.005772,Rubus_argutus:0.025281)1.000000:0.01353,Filipendula_ulmaria:0.06931)1.000000:0.006816,((((Purshia_tridentata:0.015152,Cercocarpus_ledifolius:0.012338)0.554842:0.001428,Chamaebatia_foliolosa:0.01573)0.730637:0.003057,Dryas_ajanensis_subsp._ajanensis:0.006319)1.000000:0.024606,((Aruncus_dioicus:0.026018,Spiraea_blumei:0.031786)1.000000:0.02132,((Physocarpus_opulifolius:0.013868,Neillia_sinensis:0.010046)0.999755:0.017872,((Prunus_padus:0.037371,Lyonothamnus_floribundus:0.043746)0.984562:0.002634,(((((Coleogyne_ramosissima:0.02246,Kerria_japonica:0.020587)1.000000:0.008742,Rhodotypos_scandens:0.017059)1.000000:0.00599,((Exochorda_racemosa_subsp._serratifolia:0.006462,Oemleria_cerasiformis:0.005064)0.846836:0.003373,Prinsepia_uniflora:0.018574)1.000000:0.014841)0.775411:0.004432,Sorbaria_kirilowii_var._arborea:0.044625)0.735970:0.001451,(Gillenia_trifoliata:0.016586,(Lindleya_mespiloides:0.026981,(Kageneckia_oblonga:0.033732,(Vauquelinia_australis:0.026557,(Chamaemeles_coriacea:0.023532,((((Hesperomeles_cuneata:0.002547,Phippsiomeles:0.015557)0.432917:0.013273,Crataegus_hupehensis:0.011166)0.762088:0.001276,((Malacomeles_denticulata:0.002457,Amelanchier_laevis:0.002225)1.000000:0.006092,Peraphyllum_ramosissimum:0.007237)0.998205:0.003799)0.949291:0.002889,(Dichotomanthes_tristaniicarpa:0.015388,((((((((((Karpatiosorbus_bristoliensis:0.000335,Aria_edulis:0.000028)1.000000:0.004864,((Micromeles_alnifolia:0.000663,Alniaria_alnifolia:0.001232)1.000000:0.007061,(Griffitharia_hemsleyi:0.006669,Thomsonaria_caloneura:0.031365)0.607050:0.000596)0.646404:0.001012)0.999220:0.002766,Torminalis_glaberrima:0.006015)0.999999:0.00393,((Malus_domestica:0.005315,Macromeles_tschonoskii:0.004669)1.000000:0.003499,Pyrus_communis:0.009331)0.543747:0.002135)0.976438:0.001352,(Aronia_arbutifolia:0.007589,Sorbus_aucuparia:0.003425)0.880970:0.000959)0.428566:0.00023,(Osteomeles_schweriniae:0.015001,(Pourthiaea_amphidoxa:0.00775,(Chaenomeles_speciosa:0.012251,(Pseudocydonia_sinensis:0.007781,Cydonia_oblonga:0.006948)1.000000:0.005515)0.504948:0.008304)0.845435:0.000562)0.526691:0.000252)0.482624:0.000399,(Rhaphiolepis_ferruginea:0.003729,Eriobotrya_japonica:0.000452)1.000000:0.009279)0.475066:0.000699,((Cotoneaster_frigidus:0.012316,Photinia_prunifolia:0.015522)0.672499:0.000551,Pyracantha_coccinea:0.014094)0.558990:0.000479)0.389007:0.004961,(Stranvaesia_nussia:0.00404,Weniomeles_bodinieri:0.000733)1.000000:0.008363)0.375404:0.002724,Cormus_domestica:0.026266)0.440220:0.000707)0.465065:0.000784)0.408594:0.004821)1.000000:0.005057)0.616309:0.001178)0.402887:0.002867)1.000000:0.005948)1.000000:0.016616)0.640654:0.002807)0.888602:0.001953)0.439047:0.001194)1.000000:0.003611)0.510466:0.005711)1.000000:0.067109,(Elaeagnus_angustifolia:0.223104,(Zelkova_schneideriana:0.112805,Morus_indica:0.111184)1.000000:0.059108)1.000000:0.067109);"


spname <- species.name(mptree1)
nodematrix <- read.tree.nodes(str=mptree1, name=spname)$nodes

theta_levels <- c(0.001,0.055,0.840,8.845)   # 不同程度的 ILS（从弱到强）
n_gene <- 3000                    # 每组生成3000棵基因树

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

