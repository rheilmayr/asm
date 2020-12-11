clear all
clear matrix
clear mata
set matsize 11000
set maxvar 10000
eststo clear
set more off

loc dir D:/dev/paper_repos/asm/data/
insheet using `dir'long.csv

drop if dist_amb<-300
drop if dist_aml>300
drop if year<=2001
drop if year==2017

gen y0607 = (year>2005) & (year<2008)

gen post_2003 = year>2003
gen post_2005 = year>2005
gen post_2007 = year>2007
egen long f_propid = group(propid)
egen long f_ptid = group(ptid)
egen long f_state = group(state)
egen long f_munic = group(munic)
egen prop_state = mode(state), min by(f_propid) 
egen prop_munic = mode(munic), min by(f_propid) 
egen f_prop_state = group(prop_state)
egen f_prop_munic = group(prop_munic)

sum dist_amb if dist_amb>-300 & dist_aml<300 & biome==0
gen proximity = abs(dist_amb/r(max) - 1) if dist_amb>0 & dist_aml<300
gen close = dist_amb<100 & dist_amb>0
egen gts_ever = max(gts), by(f_propid) 

label var biome "Amazon biome"
label var legal_amazon "Legal Amazon"
label var soy_suit "Suitable for soy"
label var temp "Temperature"
label var trmm "Precipitation"
label var post_2005 "Post-ASM"
label var post_2007 "Post-2007"
// label var post_2006 "Post-2006"
label var post_2003 "Post-2003"
label define biome_labels 0 "Cerrado biome" 1 "Amazon biome"
label values biome biome_labels
label define l_labels 0 "Not legal Amazon" 1 "Legal Amazon"
label values legal_amazon l_labels
label define p5_labels 0 "Pre-2005" 1 "Post-ASM"
label values post_2005 p5_labels
label define p7_labels 0 "Pre-2007" 1 "Post-2007"
label values post_2007 p7_labels

// label define p_labels 0 "Pre-2006" 1 "Post-2006"
// label values post_2006 p_labels
label define p3_labels 0 "Pre-2003" 1 "Post-2003"
label values post_2003 p3_labels
label define soy_labels 0 "Not suitable" 1 "Suitable for soy"
label values soy_suit soy_labels
label var asm_now "Amazon biome, post-2006"
label var gts_now "GTS"
label var car_now "CAR"
label define asm_labels 0 "Not ASM" 1 "Amazon biome, post-2006"
label values asm_now asm_labels
label define gts_labels 0 "Not GTS" 1 "GTS"
label values gts_now gts_labels
label define car_labels 0 "Not CAR" 1 "CAR"
label values car_now car_labels
label var proximity "Proximity"
label define close_labels 0 "Not close" 1 "Close"
label values close close_labels
label define midyear_labels 0 "not mid-years" 1 "2006-2007"
label values y0607 midyear_labels
label define prod_labels 0 "Not PRODES" 1 "PRODES"
label values prodes_mon prod_labels


///////// Table S1A: Summary stats
reg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)
estpost summarize mb2_vdefor biome legal_amazon soy_suit temp trmm roaddist urbandist set pa
esttab . using `dir'tables\ts1_summary.tex, cells("mean(fmt(a3)) sd(fmt(a3)) min(fmt(a3)) max(fmt(a3))") nomtitle nonumber fragment ///
	varlabels(mb2_vdefor "Deforested (binary)" biome "In Amazon biome (binary)" legal_amazon "In Legal Amazon (binary)" ///
	soy_suit "Suitable for soy (binary)" temp "Mean temperature (degrees C)" trmm "Annual precipitation (cm/y)" ///
	roaddist "Distance to nearest road (km)" urbandist "Distance to nearest town (km)" set "In settlement (binary)" ///
	pa "In protected area (binary)") labcol2("MapBiomas" "MMA, 2018" "MMA, 2018" ///
	"Soares-Filho et al., 2014; FAO and IIASA, 2012" "MODIS" "TRMM" "IBGE, 2010" "DNIT, 2018" ///
	"INCRA, 2018" "Funai, 2018; MMA, 2018", title(Source)) stats() replace

/// Table 1: Triple diffs
eststo clear
eststo dd_ss_inbio: areg mb2_defor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & biome==1, ///
	absorb(municcode) cluster(municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Within Amazon biome"
	
eststo dd_ss_outbio: areg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & biome==0, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"
	
eststo dd_bio_inss: areg mb2_vdefor i.biome##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & soy_suit==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Within soy-suitable"
	
eststo dd_bio_outss: areg mb2_vdefor i.biome##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & soy_suit==0, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside soy-suitable"


eststo ddd: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "All points"

	
// Export table
esttab dd_ss_inbio dd_ss_outbio dd_bio_inss dd_bio_outss ddd ///
	using `dir'tables\t1_ddd.tex, replace label se nodepvars fragment wrap ///
	keep(1.soy_suit#1.post_2005 1.biome#1.post_2005 1.biome#1.soy_suit#1.post_2005) ///
	booktabs width(0.8\hsize) alignment(c) ///
	title(Linear regression results\label{tab1}) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(sample n_points N_clust, labels("Sample" "N. points" "N. municipalities")) ///
	mtitles("Amazon DD" "Cerrado DD" "Soy-suitable DD" "Non-soy suitable DD" "Triple difference")

esttab dd_ss_inbio dd_ss_outbio dd_bio_inss dd_bio_outss ddd ///
	using `dir'tables\t1_ddd.csv, replace se ///
	keep(1.soy_suit#1.post_2005 1.biome#1.post_2005 1.biome#1.soy_suit#1.post_2005) ///
	nostar

esttab dd_ss_inbio dd_ss_outbio dd_bio_inss dd_bio_outss ddd ///
	using `dir'tables\t1_ddd_nf.csv, replace label se nodepvars wrap ///
	keep(1.soy_suit#1.post_2005 1.biome#1.post_2005 1.biome#1.soy_suit#1.post_2005) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(sample n_points N_clust, labels("Sample" "N. points" "N. municipalities")) ///
	mtitles("Amazon DD" "Cerrado DD" "Soy-suitable DD" "Non-soy suitable DD" "Triple difference")


//// Table 2: Policy interactions
gen on_car = !missing(propid)
gen on_car2 = !missing(car_year)
replace on_car2 = . if (state!="MT" & state!="PA")

scalar drop _all
program drop _all 
program define prog1, rclass
	scalar p = 2*ttail(r(df),abs(r(estimate)/r(se)))
	scalar se = string(`r(se)', "%07.5f")
	scalar t = string(`r(estimate)' / `r(se)', "%03.2f")
	scalar value = string(`r(estimate)',"%08.5f")
	if `=scalar(p)'<0.01 {
		scalar coef = "`=scalar(value)'***"
		}
	else if `=scalar(p)'<0.05 {
		scalar coef = "`=scalar(value)'**"
		}
	else if `=scalar(p)'<0.1 {
		scalar coef = "`=scalar(value)'*"
		}
	else {
		scalar coef = "`=scalar(value)'"
		}
	return local t `=scalar(t)'
	return local se `=scalar(se)'
	return local coef `=scalar(coef)'
	return local p `=scalar(p)'
end	

areg mb2_vdefor i.asm_now i.biome i.year ///
	temp trmm urbandist roaddist i.pa i.set if soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_0
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'
estadd local a_p `r(p)'

areg mb2_vdefor i.asm_now##i.car_now i.biome##i.on_car2 i.year ///
	temp trmm urbandist roaddist i.pa i.set if soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_1
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'
estadd local a_p `r(p)'
lincom 1.car_now
prog1
estadd local c `r(coef)'
estadd local c_se `r(se)'
estadd local c_p `r(p)'
lincom 1.asm_now + 1.car_now + 1.asm_now#1.car_now
prog1
estadd local ac `r(coef)'
estadd local ac_se `r(se)'
estadd local ac_p `r(p)'

areg mb2_vdefor i.asm_now##i.gts_now i.biome##i.gts_ever i.year##i.f_state ///
	temp trmm urbandist roaddist i.pa i.set if year>2001 & soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_2
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'
estadd local a_p `r(p)'
lincom 1.asm_now + 1.gts_now + 1.asm_now#1.gts_now
prog1
estadd local ag `r(coef)'
estadd local ag_se `r(se)'
estadd local ag_p `r(p)'

areg mb2_vdefor i.asm_now##i.gts_now##i.car_now i.biome##i.gts_ever##i.on_car2 i.year##i.f_state ///
	temp trmm urbandist roaddist i.pa i.set if year>2001 & soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_3
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'
estadd local a_p `r(p)'
lincom 1.car_now
prog1
estadd local c `r(coef)'
estadd local c_se `r(se)'
estadd local c_p `r(p)'
lincom 1.asm_now + 1.gts_now + 1.asm_now#1.gts_now
prog1
estadd local ag `r(coef)'
estadd local ag_se `r(se)'
estadd local ag_p `r(p)'
lincom 1.asm_now + 1.car_now + 1.asm_now#1.car_now
prog1
estadd local ac `r(coef)'
estadd local ac_se `r(se)'
estadd local ac_p `r(p)'
lincom 1.asm_now + 1.gts_now + 1.car_now + 1.asm_now#1.car_now + 1.gts_now#1.car_now + 1.asm_now#1.gts_now#1.car_now + 1.asm_now#1.gts_now
prog1
estadd local acg `r(coef)'
estadd local acg_se `r(se)'
estadd local acg_p `r(p)'


// Export table
esttab interact_0 interact_1 interact_2 interact_3 using `dir'tables\t2_complementary_policies.tex, replace label se nodepvars nomtitles fragment ///
	keep(1.asm_now) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(a a_se c c_se ac ac_se ag ag_se acg acg_se n_points N_clust, ///
	layout(@ (@) @ (@) @ (@) @ (@) @ (@) @ @) ///
	labels("ASM only" " " "CAR only" " " "ASM and CAR" " " "ASM and GTS" " " "ASM, CAR and GTS" " " ///
	"N. points" "N. municipalities"))

	
esttab interact_0 interact_1 interact_2 interact_3 using `dir'tables\t2.csv, replace p nodepvars nomtitles ///
	keep(1.asm_now) ///
	stats(a a_p c c_p ac ac_p ag ag_p acg acg_p n_points N_clust) ///
	nostar	

esttab interact_0 interact_1 interact_2 interact_3 using `dir'tables\t2_nf.csv, replace p nodepvars nomtitles ///
	keep(1.asm_now) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(a a_se c c_se ac ac_se ag ag_se acg acg_se n_points N_clust, ///
	layout(@ (@) @ (@) @ (@) @ (@) @ (@) @ @) ///
	labels("ASM only" " " "CAR only" " " "ASM and CAR" " " "ASM and GTS" " " "ASM, CAR and GTS" " " ///
	"N. points" "N. municipalities"))	
	
//// T3 - LEAKAGE
// traditional distance leakage analysis
areg mb2_vdefor c.proximity##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & biome==0 & soy_suit==1, ///
	absorb(municcode) cluster(municcode)
eststo leak_1
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"
	
areg mb2_vdefor i.post_2005##i.close i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if biome==0 & legal_amazon==1 & soy_suit==1, ///
	absorb(municcode) cluster(municcode)
eststo leak_2
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"

areg mb2_vdefor i.post_2005##i.biome##i.soy_suit i.post_2005##i.close##i.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
eststo leak_3
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"

// ILUC - Deforestation for non-soy uses. 
areg mb2_vdefor i.post_2005##i.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & a_soy_2017==0, ///
	absorb(municcode) cluster(municcode)
eststo leak_4
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.post_2005##i.biome i.post_2005##i.close i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & a_soy_2017==0, ///
	absorb(municcode) cluster(municcode)
eststo leak_5
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

// Outside of gts monitored areas
areg mb2_vdefor i.post_2005##i.biome##i.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & gts_ever==0, ///
	absorb(municcode) cluster(municcode)
eststo leak_6
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Soy-suitable, outside GTS"

// On-farm leakage - Soy-suitable properties
egen soy_prop = max(soy_suit), by(propid)
replace soy_prop = 0 if missing(propid)
areg mb2_vdefor i.post_2005##i.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & soy_prop==1, ///
	absorb(municcode) cluster(municcode)
eststo leak_7
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.post_2005##i.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & a_soy_2017==0 & soy_prop==1, ///
	absorb(municcode) cluster(municcode)
eststo leak_8
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

// On-farm leakage - original soy properties
areg mb2_vdefor i.post_2005##i.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & soy_prop_00==1, ///
	absorb(municcode) cluster(municcode)
eststo leak_9
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.post_2005##i.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & a_soy_2017==0 & soy_prop_00==1, ///
	absorb(municcode) cluster(municcode)
eststo leak_10
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

esttab leak_3 leak_2 leak_1 leak_4 leak_5 leak_6 leak_7 leak_8 leak_9 leak_10 ///
	using `dir'tables\t3_leak.tex, replace label se nodepvars fragment nonumbers wrap ///
	keep(1.post_2005#c.proximity 1.post_2005#1.close 1.post_2005#1.close#1.soy_suit 1.post_2005#1.biome 1.post_2005#1.biome#1.soy_suit) ///
	order(1.post_2005#1.biome 1.post_2005#1.biome#1.soy_suit 1.post_2005#1.close 1.post_2005#1.close#1.soy_suit 1.post_2005#c.proximity) ///
	booktabs width(0.8\hsize) alignment(c) ///
	title(Linear regression results\label{tab1}) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(n_points N_clust, labels("N. points" "N. municipalities")) ///
	mtitles("\shortstack{\\ \\ (1)}" "(2)" "(3)" "(4)" "(5)" "(6)" "\shortstack{All\\points\\(7)}" "\shortstack{Not converted\\to soy\\(8)}" "\shortstack{All\\points\\(9)}" "\shortstack{Not converted\\to soy\\(10)}") ///
	mgroups("All points" "Soy-suitable" "\shortstack{Not converted\\to soy}" "\shortstack{Not GTS-\\monitored}" "\shortstack{Properties with\\soy-suitable land}" "\shortstack{Properties growing\\soy in 2000}", ///
	pattern(1 1 0 1 0 1 1 0 1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
	
esttab leak_3 leak_2 leak_1 leak_4 leak_5 leak_6 leak_7 leak_8 leak_9 leak_10 ///
	using `dir'tables\t3_leak_nf.csv, replace label se nodepvars nonumbers wrap ///
	keep(1.post_2005#c.proximity 1.post_2005#1.close 1.post_2005#1.close#1.soy_suit 1.post_2005#1.biome 1.post_2005#1.biome#1.soy_suit) ///
	order(1.post_2005#1.biome 1.post_2005#1.biome#1.soy_suit 1.post_2005#1.close 1.post_2005#1.close#1.soy_suit 1.post_2005#c.proximity) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(n_points N_clust, labels("N. points" "N. municipalities"))
	
		
////// Figure 2C: Time effects
eststo time: areg mb2_vdefor b2005.year##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa  if legal_amazon==1 & dist_amb>-300 & dist_aml<300, absorb(municcode) cluster(municcode)	
esttab time ///
	using `dir'figures\f2_time_plot.csv, replace plain ///
	b(a3) ci(a3)
	
eststo time2: areg mb2_vdefor b2005.year##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & dist_amb>-100 & dist_amb<100, ///
	absorb(municcode) cluster(municcode)
esttab time2 ///
	using `dir'figures\f2_time_plot_100.csv, replace plain ///
	b(a3) ci(a3)
	

	
////// Figure SIB: Sensitivity to monitoring
// Within legal amazon, heterogeneous effect
eststo pol_rob_1: areg mb2_vdefor i1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if legal_amazon==1,  ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

//
eststo pol_rob_2: areg mb2_vdefor i1.y0607##i1.biome##i.soy_suit i1.post_2007##i1.biome##i.soy_suit ///
	i.post_2007##i1.prodes_mon##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

//
eststo pol_rob_3: areg mb2_vdefor 1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if ble_muni==1 & legal_amazon==1,  ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

//
eststo pol_rob_4: areg mb2_vdefor 1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if bl_muni==1 & legal_amazon==1,  ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

// All points, heterogeneous effect
eststo pol_rob_5: areg mb2_vdefor i1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa,  ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

// Legal amazon interaction, heterogeneous effect
eststo pol_rob_6: areg mb2_vdefor i1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit ///
	i1.post_2003##i1.legal_amazon##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa,  ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

// 
eststo pol_rob_7: areg mb2_vdefor i1.y0607##i1.biome##i.soy_suit i1.post_2007##i1.biome##i.soy_suit ///
	i.post_2007##i1.prodes_mon##i1.soy_suit i1.post_2003##i1.legal_amazon##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa,  ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
	
esttab pol_rob_1 pol_rob_2 pol_rob_3 pol_rob_4 pol_rob_5 pol_rob_6 pol_rob_7 ///
	using `dir'tables\ts2b_ppcdam_interact.tex, replace label se nodepvars fragment wrap nomtitles ///
	keep(1.y0607#1.biome 1.y0607#1.biome#1.soy_suit 1.post_2007#1.biome 1.post_2007#1.biome#1.soy_suit ///
		 1.post_2007#1.prodes_mon 1.post_2007#1.prodes_mon#1.soy_suit ///
		 1.post_2003#1.legal_amazon 1.post_2003#1.legal_amazon#1.soy_suit) ///
	order(1.y0607#1.biome 1.y0607#1.biome#1.soy_suit 1.post_2007#1.biome 1.post_2007#1.biome#1.soy_suit ///
		 1.post_2007#1.prodes_mon 1.post_2007#1.prodes_mon#1.soy_suit ///
		 1.post_2003#1.legal_amazon 1.post_2003#1.legal_amazon#1.soy_suit) ///	
	booktabs width(0.8\hsize) alignment(c) ///
	stats(n_points N_clust, labels("N. points" "N. municipalities")) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	mgroups("Legal Amazon" "\shortstack{Priority list\\eligible municipalities}" "\shortstack{Priority list\\municipalities}" ///
	"\shortstack{All points}", pattern(1 0 1 1 1 0 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
	
esttab pol_rob_1 pol_rob_2 pol_rob_3 pol_rob_4 pol_rob_5 pol_rob_6 pol_rob_7 ///
	using `dir'tables\tsb.csv, replace se nodepvars nomtitles ///
	keep(1.y0607#1.biome 1.y0607#1.biome#1.soy_suit 1.post_2007#1.biome 1.post_2007#1.biome#1.soy_suit ///
		 1.post_2007#1.prodes_mon 1.post_2007#1.prodes_mon#1.soy_suit ///
		 1.post_2003#1.legal_amazon 1.post_2003#1.legal_amazon#1.soy_suit) ///
	order(1.y0607#1.biome 1.y0607#1.biome#1.soy_suit 1.post_2007#1.biome 1.post_2007#1.biome#1.soy_suit ///
		 1.post_2007#1.prodes_mon 1.post_2007#1.prodes_mon#1.soy_suit ///
		 1.post_2003#1.legal_amazon 1.post_2003#1.legal_amazon#1.soy_suit) ///	
	nostar

/// Table SID: Robustness to alternate suitability definitions
gen soy_suit2 = soy_suit
eststo prim: areg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez ">Medium"


drop soy_suit2
gen soy_suit2 = suit>0
eststo soares_cat: areg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez "Not included"


drop soy_suit2
gen soy_suit2 = gaezsuit > 40
eststo gaez_cat: areg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares "Not included"
estadd loc gaez ">Medium"

drop soy_suit2
gen soy_suit2 = (suit>0) & (gaezsuit>25)
eststo prim_mod: areg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez ">Moderate"

gen soy_suit3 = gaezsuit
eststo gaez_cont: areg mb2_vdefor i.biome##c.soy_suit3##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares "Not included"
estadd loc gaez "Continuous"

drop soy_suit3
gen soy_suit3 = (suit>0) * gaezsuit
eststo both_cont: areg mb2_vdefor i.biome##c.soy_suit3##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez "Continuous"

label var soy_suit2 "Suitable for soy"
label values soy_suit2 soy_labels
label var soy_suit3 "Continuous suitability"

esttab prim soares_cat gaez_cat prim_mod gaez_cont both_cont ///
	using `dir'tables\ts5_soy_robustness.tex, replace se label wrap nomtitles fragment ///
	keep(1.biome#1.post_2005#1.soy_suit2 1.biome#1.post_2005#c.soy_suit3) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) modelwidth(8) ///
	stats(soares gaez n_points N_clust, labels("Soares-Filho et al data treatment" "GAEZ data treatment" "N. points" "N. municipalities"))

drop soy_suit2 soy_suit3

// Table SIE: Spatial and temporal subsets
eststo r1: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & dist_amb<300 & dist_amb>-300, ///
	absorb(municcode) cluster(municcode)
count if (year==2005 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc geog "Full sample"
estadd loc time "2002-2016"

eststo r2: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 ///
	& year>2003 & year<2008, ///
	absorb(municcode) cluster(municcode)
count if (year==2005 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc geog "Full sample"
estadd loc time "2004-2007"

eststo r3: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 ///
	& dist_amb>-50 & dist_amb<50, ///
	absorb(municcode) cluster(municcode)
count if (year==2005 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc geog "<50km"
estadd loc time "2002-2016"

eststo r4: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 ///
	& dist_amb>-100 & dist_amb<100, ///
	absorb(municcode) cluster(municcode)
count if (year==2005 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc geog "<100km"
estadd loc time "2002-2016"

eststo r5: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 ///
	& dist_amb>-200 & dist_amb<200, ///
	absorb(municcode) cluster(municcode)
count if (year==2005 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc geog "<200km"
estadd loc time "2002-2016"
	
eststo r6: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 ///
	& year>2003 & year<2008 & dist_amb>-100 & dist_amb<100, ///
	absorb(municcode) cluster(municcode)
count if (year==2005 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc geog "<100km"
estadd loc time "2004-2007"
	
esttab r1 r2 r3 r4 r5 r6 ///
	using `dir'tables\ts8_subsample_robustness.tex, replace se label wrap fragment ///
	keep(1.biome#1.soy_suit#1.post_2005) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) modelwidth(8) ///
	nomtitles ///
	stats(geog time n_points N_clust, labels("Spatial subsample (distance to biome boundary)" "Years included" "N. points" "N. municipalities"))
	
	
// Table SIF: Heterogeneity in cross-biome leakage (compare to Moffette and Gibbs results)
areg mb2_vdefor i.post_2005##i.close i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if biome==0 & legal_amazon==1 & soy_suit==1, ///
	absorb(municcode) cluster(municcode)
eststo state_leak_all	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.post_2005##i.close i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if state=="MT" & biome==0 & legal_amazon==1 & soy_suit==1, ///
	absorb(municcode) cluster(municcode)
eststo state_leak_mt
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.post_2005##i.close i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if state!="MT" & biome==0 & legal_amazon==1 & soy_suit==1, ///
	absorb(municcode) cluster(municcode)
eststo state_leak_notmt	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)	

esttab state_leak_all state_leak_mt state_leak_notmt ///
	using `dir'tables\ts6_leakage_by_state.tex, replace se label wrap fragment ///
	keep(1.post_2005 1.close 1.post_2005#1.close) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) modelwidth(8) ///
	mtitles("All points" "Inside Mato Grosso" "Outside Matto Grosso") ///
	stats(n_points N_clust, labels("N. points" "N. municipalities"))



// Table SIG: Tests of cross-biome leakage in different regions
gen matopiba = inlist(state, "MA", "TO", "PI", "BA")
areg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if biome==0, ///
	absorb(municcode) cluster(municcode)
eststo cerrado_all
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if biome==0 & legal_amazon==0, ///
	absorb(municcode) cluster(municcode)
eststo cerrado_nleg
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if biome==0 & legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
eststo cerrado_leg
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
	
areg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if biome==0 & matopiba==1, ///
	absorb(municcode) cluster(municcode)
eststo cerrado_matopiba
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

areg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if biome==0 & close==1, ///
	absorb(municcode) cluster(municcode)
eststo cerrado_close
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
	
esttab cerrado_all cerrado_leg cerrado_nleg cerrado_matopiba cerrado_close ///
	using `dir'tables\ts7_leakage_cerrado.csv, replace label wrap ///
	keep(1.post_2005 1.soy_suit 1.soy_suit#1.post_2005) ///
	star(* 0.1 ** 0.05 *** 0.01) modelwidth(8) ///
	mtitles("Entire Cerrado" "\shortstack{Cerrado inside\\Legal Amazon}" "\shortstack{Cerrado outside\\Legal Amazon}" "Cerrado in Matopiba" "\shortstack{Cerrado within 100km\\of Amazon biome}") ///
	stats(n_points N_clust, labels("N. points" "N. municipalities"))
	
esttab cerrado_all cerrado_leg cerrado_nleg cerrado_matopiba cerrado_close ///
	using `dir'tables\ts7_leakage_cerrado.tex, se replace label wrap fragment ///
	keep(1.post_2005 1.soy_suit 1.soy_suit#1.post_2005) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) modelwidth(8) ///
	mtitles("Entire Cerrado" "\shortstack{Cerrado inside\\Legal Amazon}" "\shortstack{Cerrado outside\\Legal Amazon}" "Cerrado in Matopiba" "\shortstack{Cerrado within 100km\\of Amazon biome}") ///
	stats(n_points N_clust, labels("N. points" "N. municipalities"))
	

/////// Table SIC - Robustness checks
/// Primary specification	
eststo primary: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "No"
estadd loc cov "Yes"
estadd loc state_fe "Yes"
estadd loc munic_fe "No"
estadd loc biome_fe "No"
estadd loc soy_fe "No"
estadd loc prop_fe "No"
estadd loc funcform "Linear"

/// Different forest definition
eststo just_mbfor: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & mb2_for_2000==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "Yes"
estadd loc cov "Yes"
estadd loc state_fe "Yes"
estadd loc munic_fe "No"
estadd loc biome_fe "No"
estadd loc soy_fe "No"
estadd loc prop_fe "No"
estadd loc funcform "Linear"

/// No extra covariates
eststo nocov: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "No"
estadd loc cov "No"
estadd loc state_fe "Yes"
estadd loc munic_fe "No"
estadd loc biome_fe "No"
estadd loc soy_fe "No"
estadd loc prop_fe "No"
estadd loc funcform "Linear"
	
/// Municipality-year fixed effects
gen y_munic = municcode * 10000 + year
eststo muniyear: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(y_munic) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "No"
estadd loc cov "Yes"
estadd loc state_fe "No"
estadd loc munic_fe "Yes"
estadd loc biome_fe "No"
estadd loc soy_fe "No"
estadd loc prop_fe "No"
estadd loc funcform "Linear"

/// Biome-year and soy-year time trends
eststo soybioyear: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 c.year##i.soy_suit c.year##i.biome ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	absorb(municcode) cluster(municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "No"
estadd loc cov "Yes"
estadd loc state_fe "Yes"
estadd loc munic_fe "No"
estadd loc biome_fe "Yes"
estadd loc soy_fe "Yes"
estadd loc prop_fe "No"
estadd loc funcform "Linear"
	
/// Property-level fixed effects
eststo propfe: areg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_prop_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & !missing(propid), ///
	absorb(f_propid) vce(cluster f_prop_munic)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "No"
estadd loc cov "Yes"
estadd loc state_fe "Yes"
estadd loc munic_fe "No"
estadd loc biome_fe "No"
estadd loc soy_fe "No"
estadd loc prop_fe "Yes"
estadd loc funcform "Linear"

/// Logistic regression - placeholder for BIFE bias-corrected regression with proper fixed effects from R
eststo logitreg: logit mb2_vdefor i.biome##i.soy_suit##i.post_2005 ///
	if _n <10000 & legal_amazon==1, ///
	vce(cluster municcode) or
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "No"
estadd loc cov "Yes"
estadd loc state_fe "Yes"
estadd loc munic_fe "No"
estadd loc biome_fe "No"
estadd loc soy_fe "No"
estadd loc prop_fe "No"
estadd loc funcform "Logistic"

/// Cox proportional hazards
preserve

keep if year==2002 | year==2006
drop if mb2_y_defor<2002

gen survival = (mb2_y_defor - year) + 1
replace survival = 4 if year==2002 & missing(survival)
replace survival = 4 if year==2002 & mb2_y_defor>=2006
replace survival = 11 if year==2006 & missing(survival)
replace survival = . if survival <0
drop if missing(survival)

gen defor = !missing(mb2_y_defor)
replace defor = 0 if year==2002 & mb2_y_defor>=2006
replace defor = 0 if year==2006 & mb2_y_defor<2006

stset survival, failure(defor)
stcox i.biome##i.soy_suit##i.post_2005 temp trmm roaddist urbandist i.pa i.set i.f_state##i.year ///
	if legal_amazon==1, vce(cluster municcode)
eststo cox
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc restrict_for "No"
estadd loc cov "Yes"
estadd loc state_fe "Yes"
estadd loc munic_fe "No"
estadd loc biome_fe "No"
estadd loc soy_fe "No"
estadd loc prop_fe "No"
estadd loc funcform "Cox"


esttab primary just_mbfor nocov muniyear soybioyear propfe logitreg cox ///
	using `dir'tables\ts4_robustness_wcox.tex, replace label se nodepvars fragment ///
	keep(1.biome#1.soy_suit#1.post_2005) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(restrict_for cov state_fe munic_fe biome_fe soy_fe prop_fe funcform n_points N_clust, ///
	labels("Restrict to MapBiomas forests" "Includes covariates" "State by year FE" ///
	"Municipality by year FE" "Biome time trend" "Soy-suitable time trend" "Property FE" "Functional form" "N. points" "N. municipalities")) ///
	nomtitles eform(0 0 0 0 0 0 1 1) eqlabels("" "" "" "" "" "" "" "")

