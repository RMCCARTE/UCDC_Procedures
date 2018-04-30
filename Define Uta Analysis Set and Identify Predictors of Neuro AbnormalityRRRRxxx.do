/***************************************************************************/
/*************  Prepare Analysis File for Uta's Paper   ********************/
/*************   Select Pts with 2+ Neuro Assessments   ********************/
/***************************************************************************/

**************  Define the input and output directories ********************
global proc    "Z:\UCDC\Procedures\LS/"
global inpath  "Z:\UCDC\Download_Jan_2018\LS\sasv9/"
global inpath2  "Z:\UCDC\Download_Jan_2018\LS\csv/"
global working "Z:\UCDC\New_Analysis_Files_Jan_2018\Working_Files/"
global anly_stata "Z:\UCDC\New_Analysis_Files_Jan_2018/Stata/"
global anly_sas "Z:\UCDC\New_Analysis_Files_Jan_2018/SAS/"
global results "Z:\UCDC\Projects\Uta_NH4_GLN_Effect_on_Neuro/"


************* Flag records to include in Uta's Analysis File
use ${anly_stata}ls_neuro_stdscrs.dta, clear
joinby participant_id using ${anly_stata}enrl_clin_socdem_ls_pcori_dx_method, unm(master) _merge(dem_clin) 
joinby participant_id using ${anly_stata}liver_transplant_date, unm(master) _merge(lt_mrg) 
keep if protocol==5101
drop if inlist(participant_id, 101983, 104351, 104125, 102129, 105438)
replace liver_trans_date=lt_dte if liver_trans_date==.
*drop seq last_rec lt_dte
sort participant_id visit_dte
gen include=1  if visit_dte < liver_trans_date   | visit_dte==.                             /**** only include records that predate liver transplant ****/
by participant_id: gen rec_n=_n  if include==1                                              /**** sequentially number all pre LT records per person  ****/
by participant_id: gen seq=_n if neuro_age~=. & include==1                                  /**** sequentially number records with neuropsych tests  ****/
tab2 rec_n seq if rec_n==1, miss                                                            /**** generate table 1 of all records vs those w/ neuro  ****/
recode ucd_dx (1/3=1 Proximal) (4/8=3 Distal), gen(prox_dist)
replace prox_dist=2 if sympt_last==0 & ucd_dx==3
label def prox_distf 1 Proximal 2 Asympt_OTC 3 Distal, replace
label values prox_dist prox_distf
label var prox_dist "Diagnostic Status"
fre prox_dist

******* Table 1 - Compare Group with Two+ Neuropsych assessments to Full Group (Needs editing)
sort participant_id src visit_age
by participant_id src: gen last_src=1 if _n==_N
tab src if last_src==1
tabout qol_agegp2plus agegp_bl ucd_dx neonatal_onset_new sex pt_race_gp pt_ethnicity fam_lang ever_sp_ed src using "Q:\Projects_Major\RDCRN_UCDC\UCDC\UCDC QOL\QOL Paper\table1br.xls" if last_src==1, mi c(freq col) clab(N Col_%) format(0 1) replace noff(3) ptotal(single) layout(cb) stats(chi2)
sort participant_id visit_age
by participant_id: gen first_rec=1 if _n==1
tabout /*qol_agegp2plus*/ agegp_bl ucd_dx neonatal_onset_new sex pt_race_gp pt_ethnicity fam_lang ever_sp_ed src using "Q:\Projects_Major\RDCRN_UCDC\UCDC\UCDC QOL\QOL Paper\table1cr.xls" if first_rec==1, mi c(freq col) clab(N Col_%) format(0 1) replace noff(3) ptotal(single) layout(cb) stats(chi2)
 
save "Q:\Projects_Major\RDCRN_UCDC\UCDC\UCDC QOL\QOL Paper\QOL Paper Dataset.dta", replace
 
**********   Perform Analysis to Estimate QOL Scores by Source
 
use "Q:\Projects_Major\RDCRN_UCDC\UCDC\UCDC QOL\QOL Paper\QOL Paper Dataset.dta", clear
pause on
***** Child Self-Assessment
*tabstat qol_age physical_score emotion_score social_score school_score psychosocial_health physical_health qol_total if last_src==1 & src==1, by(neonatal_onset_new) stat(n mean sd median min max) long col(stat)
foreach y in physical_score emotion_score social_score school_score psychosocial_health physical_health qol_total {
qreg `y'  i.neonatal if last_src==1 & src==1
margins neonatal
}

by participant_id: gen first_rec=1 if seq==1 & include==1                                   /**** define the first pre-LT neuropsych test record     ****/
by participant_id: gen begin_dte=visit_dte[1] if include==1                                 /**** define start date of the first neuropsych record   ****/
format begin_dte %td
by participant_id: gen last_rec=1 if neuro_age~=. & seq==_N & include==1                    /**** define the last pre-lt neuropsych test record      ****/
gen uta_set = 1 if ucd_dx~=9 & protocol==5101 & (seq==1 | last_rec==1) & include==1
replace uta_set=. if last_rec==1 & first_rec==1 
*by participant_id: gen begin_dte=visit_dte[1] if uta_set==1 & last_rec==1

tab2 ucd_dx seq if last_rec==1 & uta_set==1
tab2 enrl_status seq if last_rec==1 & uta_set==1, miss
order first_rec last_rec rec_n seq neuro_age, after (participant_id)
l participant_id begin_dte visit_dte liver_trans_date if last_rec==1 & uta_set==1 & liver_trans_status==1

global glob_func "bl_fscale_iq_std bl_abas_gac_prt_sc"
global lang      "bl_verb_iq_std bl_perf_iq_std" 
global visual    "bl_beery_visual_std bl_beery_visualpercep_std"
global motor     "bl_groove_peg_d_std bl_groove_peg_n_std bl_grip_str_d_std bl_grip_str_n_std"
global att_exec  "bl_cbcl_att_prb_prt_std bl_teach_std bl_teachwalkdontwalk_std bl_brief_shft_prt_std bl_tldttlmovestd"
global mem       "bl_cvlt_lista_trial5_std bl_cvlt_shrt_dly_recall_std bl_cvlt_long_dly_recall_std bl_brief_wkmem_prt_std"
global behav     "bl_cbcl_int_prb_prt_std bl_cbcl_ext_prb_prt_std"


***************************************************************************************************************
******* Create Baseline (first neuro assessment in interval) and difference from Baseline to last assessment
***************************************************************************************************************
foreach y in bayleycognitivecomposite fscale_iq_std abas_gac_prt_sc perf_iq_std bayleylangcompcomposite verb_iq_std beery_visual_std beery_visualpercep_std bayleymotorcompcomposite groove_peg_d_std groove_peg_n_std grip_str_d_std grip_str_n_std  /// 
     cbcl_att_prb_prt_std teach_std teachwalkdontwalk_std brief_shft_prt_std tldttlmovestd cvlt_lista_trial5_std cvlt_shrt_dly_recall_std cvlt_long_dly_recall_std brief_wkmem_prt_std /// 
	  cbcl_int_prb_prt_std cbcl_ext_prb_prt_std {
gen b_`y' = `y' if seq==1 & uta_set==1
by participant_id: egen bl_`y'=mean(b_`y') 
gen dif_`y' = `y' - bl_`y'
drop b_`y'
replace bl_`y'= . if dif_`y'==.
*display ""
*display "Last vs First Difference in `y'"
*signrank dif_`y'=0 if uta_set==1 & last_rec==1
}

label def neonatef 0 "Later" 1 "Neonatal", replace 
label values neonatal_onset neonatef
label var neonatal_onset "Timing of Onset"  

****************************************************************************************************************
******* First vs Last Assessment Differences Overall - 
*******  Table 2 - Neuropsych baseline and difference by Proximal vs Distal UCD 
****************************************************************************************************************

***** For those aged < 3 yrs
foreach y in bayleycognitivecomposite bayleylangcompcomposite bayleymotorcompcomposite {
*display ""
*display "Baseline `y'"
*display ""
*anova bl_`y' i.prox_dist if last_rec==1
*bys neonatal_onset: summ dif_`y' if last_rec==1

display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' if last_rec==1, vce(robust)
qreg bl_`y' i.prox_dist if last_rec==1, vce(robust)
margins prox_dis
contrast prox_dist, asobserved
*marginsplot, recast(bar) name(bar_`y', replace) title ("Change in `y' from First to Last Assessment") ytitle(Change First to Last)
}

***** For those 3 yrs and older
foreach y in fscale_iq_std abas_gac_prt_sc  verb_iq_std perf_iq_std beery_visual_std beery_visualpercep_std groove_peg_d_std groove_peg_n_std grip_str_d_std grip_str_n_std  /// 
     cbcl_att_prb_prt_std teach_std teachwalkdontwalk_std brief_shft_prt_std tldttlmovestd cvlt_lista_trial5_std cvlt_shrt_dly_recall_std cvlt_long_dly_recall_std brief_wkmem_prt_std /// 
	  cbcl_int_prb_prt_std cbcl_ext_prb_prt_std {
*display ""
*display "Baseline `y'"
*display ""
*anova bl_`y' i.prox_dist if last_rec==1
*bys neonatal_onset: summ dif_`y' if last_rec==1

display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' if last_rec==1, vce(robust)
qreg bl_`y' i.prox_dist if last_rec==1, vce(robust)
margins prox_dist
contrast prox_dist, asobserved
*marginsplot, recast(bar) name(bar_`y', replace) title ("Change in `y' from First to Last Assessment") ytitle(Change First to Last)
}

tabstat bl_bayleycognitivecomposite bl_bayleylangcompcomposite bl_bayleymotorcompcomposite if last_rec==1 & uta_set==1 & form==1,  by(prox_dist) stat(n mean sd median min max) long col(stat)
tabstat $glob_func $lang $visual $motor $att_exec $mem $behav if last_rec==1 & uta_set==1 & inlist(formn,2, 3, 4, 5),  by(prox_dist) stat(n mean sd median min max) long col(stat)


*************************************************************************************************************
******* First vs Last Assessment Differences within Neonatal and Late Onset UCD and by whether pts had an HA event
****************************************************************************************************************

gen HA_event= 0
replace HA_event = 1 if sympt_last==2 

label def neonatef 0 "Later" 1 "Neonatal", replace
label values neonatal_onset neonatef
label var neonatal_onset "Timing of Onset" 
label def yonf 0 No 1 Yes, replace
label values HA_event yonf
label var HA_event "Experienced Hyperammonemic Event(s)"



foreach y in bayleycognitivecomposite fscale_iq_std abas_gac_prt_sc bayleylangcompcomposite verb_iq_std perf_iq_std groove_peg_d_std groove_peg_n_std grip_str_d_std grip_str_n_std beery_visual_std beery_visualpercep_std /// 
     cbcl_att_prb_prt_std /*teach_std */ teachwalkdontwalk_std brief_shft_prt_std tldttlmovestd cvlt_lista_trial5_std cvlt_shrt_dly_recall_std cvlt_long_dly_recall_std brief_wkmem_prt_std /// 
	  cbcl_int_prb_prt_std cbcl_ext_prb_prt_std {
display ""
display "Last vs First Difference in `y'"
display ""
summ `y' 
/*
qreg dif_`y' i.neonatal_onset if last_rec==1, vce(robust) iter(200)
margins neonatal
marginsplot, recast(bar) name(bar_neo_`y', replace) title ("Change in `y' from First to Last Assessment") ytitle(Change First to Last)
*/
}

foreach y in bayleycognitivecomposite fscale_iq_std abas_gac_prt_sc bayleylangcompcomposite perf_iq_std verb_iq_std groove_peg_d_std groove_peg_n_std grip_str_d_std grip_str_n_std beery_visual_std beery_visualpercep_std /// 
     cbcl_att_prb_prt_std  teach_std teachwalkdontwalk_std brief_shft_prt_std tldttlmovestd cvlt_lista_trial5_std cvlt_shrt_dly_recall_std cvlt_long_dly_recall_std brief_wkmem_prt_std /// 
	   cbcl_int_prb_prt_std cbcl_ext_prb_prt_std  {
display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' i.HA_event if last_rec==1, vce(robust) iter(200)
margins HA_event
marginsplot, recast(bar) name(bar_HA_`y', replace) title ("Change in `y' from First to Last Assessment") ytitle(Change First to Last)
}

drop if uta_set~=1 | last_rec~=1
save ${anly_stata}Uta_Analysis_DatasetR.dta, replace

****************************************************************************************************************
******* Merge record of hyperammonemic events 
****************************************************************************************************************
joinby participant_id using ${anly_stata}HA_events.dta, unm(master) _merge(ha_mrg)
fre ha_mrg
gen dte=visit_dte 
replace dte=admit_dte if admit_dte ~=.
keep if dte <= visit_dte & dte >= begin_dte & last_rec==1 & uta_set==1
fre ha_mrg
sort participant_id dte
by participant_id: gen first_ha=1 if _n==1
by participant_id: gen last_ha=1 if _n==_N
by participant_id: egen HA_cnt = count(HA_nh4_above_100)
replace HA_cnt=0 if HA_cnt==.
fre HA_cnt if last_ha==1
recode HA_cnt (0=0)(1=1)(2=2) (3/5=3 "3-5") (6/19 = 4 "6+"), gen(HA_frq)
by participant_id: egen max_nh4 = max(peak_ammonia)
by participant_id: egen max_glm = max(glut_presenting)
by participant_id: egen HA_coma = max(coma)
keep if last_ha==1
drop ha_mrg admit_dte dschg_dte HA peak_ammonia glut_presenting coma visit_type HA_nh4_above_100 dte first_ha last_ha
save ${anly_stata}Uta_Analysis_DatasetR_HA_events.dta, replace

****************************************************************************************************************
******* Merge laboratory assessments and trim to only those collected in between first and last neuro assessment
****************************************************************************************************************
joinby participant_id using ${anly_stata}standard_lab_data.dta, unm(master) _merge(lab_mrg)
fre lab_mrg
gen dte=visit_dte 
replace dte=lab_reference_date if lab_reference_date ~=.
keep if dte <= visit_dte & dte >= begin_dte & last_rec==1 & uta_set==1
sort participant_id dte
by participant_id: gen first_lab=1 if _n==1
by participant_id: gen last_lab=1 if _n==_N
tab uta_set if last_lab==1
gen time = (visit_dte - begin_dte)/30.4375 

****************************************************************************************************************
******* Create medians bouts of elevation to represent ammonia and glutamine levels in the interim between tests
****************************************************************************************************************
summ ammonia, detail
summ glutamine, detail
gen nh4_gt50= 1 if ammonia > 50 & ammonia ~=.
gen nh4_gt70= 1 if ammonia > 70 & ammonia ~=.
gen glm_gt800= 1 if glutamine > 800 & glutamine ~=.
gen glm_gt1000=1 if glutamine > 1000 & glutamine ~=.


by participant_id: egen ammonia_mid_hi=count(nh4_gt50)
by participant_id: egen ammonia_hi = count(nh4_gt70)
by participant_id: egen glutamine_mid_hi = count (glm_gt800)
by participant_id: egen glutamine_hi = count (glm_gt1000)

by participant_id: egen ammonia_p50= median(ammonia)
by participant_id: egen ammonia_p75= pctile(ammonia), p(75)
by participant_id: egen glutamine_p50 = median(glutamine) 
by participant_id: egen glutamine_p75 = pctile(glutamine), p(75)
summ ammonia_p50 glutamine_p50 if last_lab==1

table neonatal_onset if last_lab==1, row c (n ammonia_p50 p50 ammonia_p50 iqr ammonia_p50 p50 glutamine_p50 iqr glutamine_p50) 
table neonatal_onset if last_lab==1, row c (n ammonia_p75 p50 ammonia_p75 iqr ammonia_p75 p50 glutamine_p75 iqr glutamine_p75) 
table HA_event if last_lab==1, row c (n ammonia_p50 p50 ammonia_p50 iqr ammonia_p50 p50 glutamine_p50 iqr glutamine_p50) 
table HA_event if last_lab==1, row c (n ammonia_p75 p50 ammonia_p75 iqr ammonia_p75 p50 glutamine_p75 iqr glutamine_p75) 
keep if uta_set==1 & last_lab==1

save ${anly_stata}Uta_Analysis_Dataset_with_HA_events_labs.dta, replace
****************************************************************************************************************
******* Evaluate Change in Ammonia and Glutamine over the Interval between First and Last Neuropsych Eval
****************************************************************************************************************

qreg ammonia c.time##c.time, vce(robust)
margins, at(time=(0(20)140))
marginsplot, addplot(scatter ammonia time, xlabel(0(20)140) ylabel(0(50)300)) title(Ammonia Levels in the Interim between Neuropsychological Tests) subtitle(with Median Regression Curve +/- 95% CI) ///
 ytitle(Ammonia Level) xtitle(Time Interval(mos) between Neuropsychological Tests) name(All_NH4_Int, replace)

qreg ammonia_p50 c.time if last_lab==1, vce(robust)
margins, at(time=(0(20)140))
marginsplot, addplot(scatter ammonia_p50 time, xlabel(0(20)140) ylabel(0(50)300)) title(Median Ammonia Levels in the Interim between Neuropsychological Tests) subtitle(with Median Regression Curve +/- 95% CI) ///
 ytitle(Ammonia Level) xtitle(Time Interval(mos) between Neuropsychological Tests) name(P50_NH4_Int, replace)
 
qreg ammonia_p75 c.time if last_lab==1, vce(robust)
margins, at(time=(0(20)140))
marginsplot, addplot(scatter ammonia_p75 time, xlabel(0(20)140) ylabel(0(50)300)) title(Seventy-fifth Percentile Ammonia Levels in the Interim between Neuropsychological Tests) subtitle(with Median Regression Curve +/- 95% CI) ///
 ytitle(Ammonia Level) xtitle(Time Interval(mos) between Neuropsychological Tests) name(P75_NH4_Int, replace)
 
qreg glutamine c.time##c.time, vce(robust)
margins, at(time=(0(20)140))
marginsplot, addplot(scatter glutamine time, xlabel(0(20)140) ylabel(0(500)2000)) title(Glutamine Levels in the Interim between Neuropsychological Tests) subtitle(with Median Regression Curve +/- 95% CI) ///
 ytitle(Glutamine Level) xtitle(Time Interval(mos) between Neuropsychological Tests) name(All_Glut_Int, replace)
 
qreg glutamine_p50 c.time##c.time if last_lab==1, vce(robust)
margins, at(time=(0(20)140))
marginsplot, addplot(scatter glutamine_p50 time, xlabel(0(20)140) ylabel(0(500)2000)) title(Median Glutamine Levels in the Interim between Neuropsychological Tests) subtitle(with Median Regression Curve +/- 95% CI) ///
 ytitle(Glutamine Level) xtitle(Time Interval(mos) between Neuropsychological Tests) name(P50_Glut_Int, replace)
 
qreg glutamine_p75 c.time##c.time if last_lab==1, vce(robust)
margins, at(time=(0(20)140))
marginsplot, addplot(scatter glutamine_p75 time, xlabel(0(20)140) ylabel(0(500)2000)) title(Seventy-fifth Percentile Glutamine Levels in the Interim between Neuropsychological Tests) subtitle(with Median Regression Curve +/- 95% CI) ///
 ytitle(Glutamine Level) xtitle(Time Interval(mos) between Neuropsychological Tests) name(P75_Glut_Int, replace)

qreg ammonia_p50 c.glutamine_p50 c.time if last_lab==1, vce(robust)
margins, at(glutamine_p50=(300(300)1500))
marginsplot, addplot(scatter ammonia_p50 glutamine_p50, xlabel(0(300)1500) ylabel(0(50)200)) title(Ammonia vs. Glutamine Levels in the Interim between Neuropsychological Tests) subtitle(with Median Rgegression Curve +/- 95% CI) ///
 ytitle(Ammonia Level) xtitle(Glutamine Level) name(P50_NH4_Glut_Int, replace)
 
 
 
*****************************************************************************************************************
************** Analysis of relationship between metabolite levels and change in 
**************  Bayley Cog, Bayley Lang and Grip Str and CVLT Std scores
*****************************************************************************************************************


******* Analysis by presence vs absence of HA events in the interim between two neuropsych assessments

foreach y in /* bayleycognitivecomposite bayleylangcompcomposite */ beery_visual_std grip_str_d_std grip_str_n_std cvlt_lista_trial5_std {
display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' /* c.glutamine_hi##c.ammonia_hi */ i.HA_b i.sex  c.bl_`y'  c.time  c.neuro_age, vce(robust) 
margins HA_b 
marginsplot, recast(bar) title ( Change in `y' by Presence of Interim Hyperammonia)  ytitle(Change in `y') xtitle(Hyperammonemia Status) name(`y'_ha, replace)
}


******* Analysis by presence vs absence of HA events as well as Ammonia and Glutamine excusions in the interim between two neuropsych assessments
******* Beery Analysis
foreach y in /* bayleycognitivecomposite bayleylangcompcomposite */ beery_visual_std /*grip_str_d_std grip_str_n_std cvlt_lista_trial5_std */ {
display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' i.glm_hi2 i.nh4_hi2##i.HA_b  i.sex  c.bl_`y'   c.time   c.neuro_age, vce(robust) 
margins glm_hi2 
marginsplot, recast(bar) xdim(glm_hi2 ) title ( "Change in `y' by Presence of Bouts of High Glutamine, Ammonia")  ytitle(Change in `y') xtitle( "Bouts of High Glutamine, Ammonia Status") name(`y'_nh4_glm, replace)
margins nh4_hi2#HA_b
marginsplot, recast(bar) xdim(HA_b nh4_hi2) title ("Change in `y' by Presence of Interim Hyperammonia, Bouts of High Ammonia")  ytitle(Change in `y') xtitle("Hyperammonemia, Bouts of High Ammonia Status") name(`y'_ha_nh4, replace)
}

****** Analysis by presence vs absence of HA events as well as Ammonia and Glutamine excusions in the interim between two neuropsych assessments
******* Grip Strength Analysis
foreach y in /* bayleycognitivecomposite bayleylangcompcomposite beery_visual_std */ grip_str_d_std grip_str_n_std /* cvlt_lista_trial5_std */ {
display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' i.glm_hi2 i.nh4_hi2  i.HA_b  c.bl_`y'  i.sex c.time   c.neuro_age, vce(robust) 
margins HA_b 
marginsplot, recast(bar) title ( Change in `y' by Presence of Interim Hyperammonia)  ytitle(Change in `y') xtitle(Hyperammonemia Status) name(`y'_ha_nh4_glm, replace)
margins  glm_hi2 
marginsplot, recast(bar) title ( "Change in `y' by Presence of Bouts of High Glutamine, Ammonia")  ytitle(Change in `y') xtitle( "Bouts of High Glutamine Status") name(`y'_HE_nh4, replace)
margins nh4_hi2
marginsplot, recast(bar) title ("Change in `y' by Presence of Bouts of High Ammonia")  ytitle(Change in `y') xtitle("Bouts of High Ammonia Status") name(`y'_HE_glm, replace)
}


****** Analysis by presence vs absence of HA events as well as Ammonia and Glutamine excusions in the interim between two neuropsych assessments
******* CVLT Analysis
foreach y in /* bayleycognitivecomposite bayleylangcompcomposite beery_visual_std grip_str_d_std grip_str_n_std */ cvlt_lista_trial5_std {
display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' i.glm_hi2 i.HA_b i.nh4_hi2 i.HA_b  /* i.sex  */  c.bl_`y' c.time  c.neuro_age, vce(robust) 
margins HA_b 
marginsplot, recast(bar) title ( Change in `y' by Presence of Interim Hyperammonia)  ytitle(Change in `y') xtitle(Hyperammonemia Status) name(`y'_ha, replace)
margins  glm_hi2 
marginsplot, recast(bar) title ( "Change in `y' by Presence of Bouts of High Glutamine, Ammonia")  ytitle(Change in `y') xtitle( "Bouts of High Glutamine Status") name(`y'_glm, replace)
margins nh4_hi2
marginsplot, recast(bar) title ("Change in `y' by Presence of Bouts of High Ammonia")  ytitle(Change in `y') xtitle("Bouts of High Ammonia Status") name(`y'_nh4, replace)
}

foreach y in bayleycognitivecomposite bayleylangcompcomposite beery_visual_std grip_str_d_std grip_str_n_std cvlt_lista_trial5_std {
display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' c.ammonia_p50##c.bl_`y' c.time c.neuro_age i.sex  
margins, at(ammonia_p50=(0 50 100 150 200 250 300))
marginsplot, title ( Change in `y' by Interim Ammonia Level)  ytitle(Change in `y') xtitle(Interim Ammonia Level) name(`y'_nh4, replace)
} 

foreach y in bayleycognitivecomposite bayleylangcompcomposite beery_visual_std  grip_str_d_std grip_str_n_std cvlt_lista_trial5_std {
display ""
display "Last vs First Difference in `y'"
display ""
qreg dif_`y' c.glutamine_p50 c.neuro_age c.time /* c.bl_grip_str_d_std */ i.sex if last_lab==1, vce(robust)
margins, at(glutamine=(0 100 250 400 600 800 1000 1200 1400 1600 1800))
*marginsplot, title ( Change in `y' by Interim Glutamine Level)  ytitle(Change in `y') xtitle(Interim Glutamine Level) name(`y'_glut_nh4, replace)
*margins, at(ammonia=(0 50 100 150 200 250 300))
*marginsplot, title ( Change in `y' by Interim Ammonia Level)  ytitle(Change in `y') xtitle(Interim Ammonia Level) name(`y'_nh4_glut, replace)
}
keep if last_lab==1
save ${anly_stata}Uta_Analysis_Dataset_Median_Lab_Results.dta, replace

*****************************************************************************************************************
************** Analysis of relationship between metabolite levels and overall change in 
************** Bayley Cog, Beery Visual, Grip Strength and CVLT
*****************************************************************************************************************

foreach y in   /* abas_gac_prt_sc */ brief_shft_prt_std  /* brief_wkmem_prt_std cbcl_int_prb_prt_std */ cbcl_ext_prb_prt_std{
qreg dif_`y' c.glutamine##i.neonatal_onset c.bl_`y' c.time c.neuro_age i.sex if last_lab==1, vce(robust) 
margins neonatal_onset, at(glutamine=(0 100 250 400 600 800 1000 1200 1400 1600 1800))
marginsplot, title ( Change in `y' by Interim Glutamine Level)  ytitle(Change in `y') xtitle(Interim Glutamine Level) name(`y'_glm_nh4, replace)
}
*****************************************************************************************************************
************** Analysis of relationship between metabolite levels and change in ABAS, Brief, CBCL 
*****************************************************************************************************************

foreach y in   /* abas_gac_prt_sc */ brief_shft_prt_std  /* brief_wkmem_prt_std cbcl_int_prb_prt_std */ cbcl_ext_prb_prt_std{
qreg dif_`y' c.glutamine##i.neonatal_onset c.bl_`y' c.time c.neuro_age i.sex if last_lab==1, vce(robust) 
margins neonatal_onset, at(glutamine=(0 100 250 400 600 800 1000 1200 1400 1600 1800))
marginsplot, title ( Change in `y' by Interim Glutamine Level)  ytitle(Change in `y') xtitle(Interim Glutamine Level) name(`y'_glm_nh4, replace)
}

foreach y in  abas_gac_prt_sc /* brief_shft_prt_std brief_wkmem_prt_std */ cbcl_int_prb_prt_std cbcl_ext_prb_prt_std {
qreg dif_`y' c.ammonia##c.ammonia##i.sympt_last c.bl_`y' c.time c.neuro_age i.sex if last_lab==1, vce(robust) 
margins sympt_last, at(ammonia=(0 50 100 150 200 250 300))
marginsplot, title ( Change in `y' by Interim Ammonia Level)  ytitle(Change in `y') xtitle(Interim Ammonia Level) name(`y'_nh4, replace)
} 



*****************************************************************************************************
********************* Older Analyses 
****************************************************************************************************

/*
***************************************** Create Ammonia and Glutamine Categories
summ glutamine, detail
summ ammonia, detail
egen ammonia_gp = cut(ammonia), g(5) label
egen glutamine_gp = cut(glutamine), g(5) label
fre ammonia_gp glutamine_gp
gen amm_glut = ammonia * glutamine

foreach x in  glutamine{
*gen `x'_time = `x' * time 
foreach y in grip_str_d_std grip_str_n_std cvlt_lista_trial5_std {
*qreg2  `x' time if last_lab==1, cluster(participant_id) 
qreg dif_`y' `x' bl_`y' time if last_lab==1, vce(robust)

}
} 

*/
/*
qreg2 dif_cvlt_lista_trial5_std glutamine bl_cvlt_lista_trial5_std time, cluster(participant_id) 
adjust, se ci by(time) gen(p_cvlt_d e_cvlt_d)
graph twoway (lfitci p_cvlt_d glutamine, xlabel(0 100 250 400 600 800 1000 1200 1400 1600 1800) name(CVLT_Glut, replace))/*|| (lfitci p`x' visit if gp==1 & time_mos <= 84), xlabel(6 12 24 36 48 60 72) */
drop p_cvlt_d e_cvlt_d
*/



*qreg2 dif_grip_str_d_std ammonia glutamine bl_grip_str_d_std time, cluster(participant_id) 

/*
foreach y in cbcl_int_prb_prt_std cbcl_ext_prb_prt_std {
reg dif_`y' bl_`y' i.neonatal i.sympt_last if last_rec==1, robust 
margins neonatal sympt_last
reg `y' bl_`y' i.neonatal i.sympt_last if last_rec==1, robust 
margins neonatal sympt_last
}
*/
