clear all
clear matrix
clear mata
set matsize 11000
set maxvar 10000
eststo clear
set more off

loc dir D:\cloud\Dropbox\collaborations\glue-sb\soyM\analysis\8-20-19\
insheet using `dir'long.csv

// replace gts = 0 if year==2007 //// NOTE GTS VARIABLE CREATION SHOULD BE CORRECTED IN DATA GENERATION FILE
// drop gts_now
// gen gts_now = prodes_mon * gts * biome * !set * !pa
// drop soy_suit
// gen soy_suit = (suit>0 & gaezsuit>40)
// gen mb_defor = mb2_vdefor


drop if dist_amb<-300
drop if dist_aml>300
drop if year<=2001

gen y0607 = (year>2005) & (year<2008)

gen post_2003 = year>2003
// gen post_2004 = year>2004 // Change to 2003?
gen post_2005 = year>2005
// gen post_2006 = year>2006
gen post_2007 = year>2007
// gen post_2008 = year>2008
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


///////// Table S1: Summary stats
reg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)
estpost summarize mb2_vdefor biome legal_amazon soy_suit temp trmm roaddist urbandist set pa
esttab using `dir'tables\ts1_summary.tex, cells("mean(fmt(a3)) sd(fmt(a3)) min(fmt(a3)) max(fmt(a3))") nomtitle nonumber fragment ///
	varlabels(mb2_vdefor "Deforested (binary)" biome "In Amazon biome (binary)" legal_amazon "In Legal Amazon (binary)" ///
	soy_suit "Suitable for soy (binary)" temp "Mean temperature (degrees C)" trmm "Annual precipitation (cm/y)" ///
	roaddist "Distance to nearest road (km)" urbandist "Distance to nearest town (km)" set "In settlement (binary)" ///
	pa "In protected area (binary)") labcol2("Mapbiomas" "MMA, 2018" "MMA, 2018" ///
	"Soares-Filho et al., 2014; FAO and IIASA, 2012" "MODIS" "TRMM" "IBGE, 2010" "DNIT, 2018" ///
	"INCRA, 2018" "Funai, 2018; MMA, 2018", title(Source)) stats() replace

/// Table 1: Triple diffs
eststo clear
eststo dd_ss_inbio: reg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & biome==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Within Amazon biome"
	
eststo dd_ss_outbio: reg mb2_vdefor i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & biome==0, ///
	vce(cluster municcode)		
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"
	
eststo dd_bio_inss: reg mb2_vdefor i.biome##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & soy_suit==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Within soy-suitable"
	
eststo dd_bio_outss: reg mb2_vdefor i.biome##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & soy_suit==0, ///
	vce(cluster municcode)		
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside soy-suitable"

eststo ddd: reg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
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

// 	indicate(Year by state effects = *.year) ///
// 	addnotes("Includes additional covariates and interaction terms as described in SI. Standard errors clustered by municipality") ///

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
end	

reg mb2_vdefor i.asm_now i.biome i.year ///
	temp trmm urbandist roaddist i.pa i.set if soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_0
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'

reg mb2_vdefor i.asm_now##i.car_now i.biome##i.on_car2 i.year ///
	temp trmm urbandist roaddist i.pa i.set if soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_1
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'
lincom 1.car_now
prog1
estadd local c `r(coef)'
estadd local c_se `r(se)'
lincom 1.asm_now + 1.car_now + 1.asm_now#1.car_now
prog1
estadd local ac `r(coef)'
estadd local ac_se `r(se)'

reg mb2_vdefor i.asm_now##i.gts_now i.biome##i.gts_ever i.year##i.f_state ///
	temp trmm urbandist roaddist i.pa i.set if year>2001 & soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_2
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'
lincom 1.asm_now + 1.gts_now + 1.asm_now#1.gts_now
prog1
estadd local ag `r(coef)'
estadd local ag_se `r(se)'

reg mb2_vdefor i.asm_now##i.gts_now##i.car_now i.biome##i.gts_ever##i.on_car2 i.year##i.f_state ///
	temp trmm urbandist roaddist i.pa i.set if year>2001 & soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo interact_3
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_se `r(se)'
lincom 1.car_now
prog1
estadd local c `r(coef)'
estadd local c_se `r(se)'
lincom 1.asm_now + 1.gts_now + 1.asm_now#1.gts_now
prog1
estadd local ag `r(coef)'
estadd local ag_se `r(se)'
lincom 1.asm_now + 1.car_now + 1.asm_now#1.car_now
prog1
estadd local ac `r(coef)'
estadd local ac_se `r(se)'
lincom 1.asm_now + 1.gts_now + 1.car_now + 1.asm_now#1.car_now + 1.gts_now#1.car_now + 1.asm_now#1.gts_now#1.car_now + 1.asm_now#1.gts_now
prog1
estadd local acg `r(coef)'
estadd local acg_se `r(se)'


// Export table
loc dir D:\cloud\Dropbox\collaborations\glue-sb\soyM\analysis\8-20-19\
esttab interact_0 interact_1 interact_2 interact_3 using `dir'tables\t2_complementary_policies.tex, replace label se nodepvars nomtitles fragment ///
	keep(1.asm_now) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(a a_se c c_se ac ac_se ag ag_se acg acg_se n_points N_clust, ///
	layout(@ (@) @ (@) @ (@) @ (@) @ (@) @ @) ///
	labels("ASM only" " " "CAR only" " " "ASM and CAR" " " "ASM and GTS" " " "ASM, CAR and GTS" " " ///
	"N. points" "N. municipalities"))

// 	addnotes("Models estimated using forested, soy-suitable points in the states of Mato Grosso and Para. Models include additional covariates and interaction terms as described in SI. Standard errors were clustered by municipalities." ///
// 	"Terms presented in this table represent linear combinations of coefficiencts on individual and interaction terms to capture aggregate effect of multiple policies.") ///
	
//// T3 - LEAKAGE
// traditional distance leakage analysis
reg mb2_vdefor c.proximity##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & biome==0 & soy_suit==1, ///
	vce(cluster municcode)
eststo leak_1
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"
	
reg mb2_vdefor i.post_2005##i.close i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if biome==0 & legal_amazon==1 & soy_suit==1, ///
	vce(cluster municcode)	
eststo leak_2
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"

reg mb2_vdefor i.post_2005##i.biome##i.soy_suit i.post_2005##i.close##i.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1, ///
	vce(cluster municcode)
eststo leak_3
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"

// Outside of gts monitored areas
reg mb2_vdefor i.post_2005##i.biome##i.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & gts_ever==0, ///
	vce(cluster municcode)	
eststo leak_4
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Soy-suitable, outside GTS"

// ILUC - Deforestation for non-soy uses. 
reg mb2_vdefor i.post_2005##i.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & a_soy_2017==0, ///
	vce(cluster municcode)
eststo leak_5
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

reg mb2_vdefor i.post_2005##i.biome i.post_2005##i.close i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set ///
	if legal_amazon==1 & a_soy_2017==0, ///
	vce(cluster municcode)
eststo leak_6
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

loc dir D:\cloud\Dropbox\collaborations\glue-sb\soyM\analysis\8-20-19\
esttab leak_3 leak_2 leak_1 leak_4 leak_5 leak_6 ///
	using `dir'tables\t3_leak.tex, replace label se nodepvars nomtitles fragment ///
	keep(1.post_2005#c.proximity 1.post_2005#1.close 1.post_2005#1.close#1.soy_suit 1.post_2005#1.biome 1.post_2005#1.biome#1.soy_suit) ///
	order(1.post_2005#1.biome 1.post_2005#1.biome#1.soy_suit 1.post_2005#1.close 1.post_2005#1.close#1.soy_suit 1.post_2005#c.proximity) ///
	booktabs width(0.8\hsize) alignment(c) ///
	title(Linear regression results\label{tab1}) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(n_points N_clust, labels("N. points" "N. municipalities")) ///
	mgroups("All points" "Soy-suitable" "\shortstack{Not GTS-\\monitored}" "\shortstack{Not converted\\to soy}" ///
	"\shortstack{All points}", pattern(1 1 0 1 1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span}))
		
// 	indicate(Year by state effects = *.year) ///
// 	addnotes("Models estimated using forested points located within the Legal Amazon. Models include additional covariates and interaction terms as described in SI. Standard errors clustered by municipality.")
	
	
////// Figure SIB: Sensitivity to monitoring
// Within legal amazon, heterogeneous effect
eststo pol_rob_1: reg mb2_vdefor i1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if legal_amazon==1, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

//
eststo pol_rob_2: reg mb2_vdefor i1.y0607##i1.biome##i.soy_suit i1.post_2007##i1.biome##i.soy_suit ///
	i.post_2007##i1.prodes_mon##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if legal_amazon==1, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

//
eststo pol_rob_3: reg mb2_vdefor 1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if ble_muni==1 & legal_amazon==1, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

//
eststo pol_rob_4: reg mb2_vdefor 1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if bl_muni==1 & legal_amazon==1, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

// All points, heterogeneous effect
eststo pol_rob_5: reg mb2_vdefor i1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

// Legal amazon interaction, heterogeneous effect
eststo pol_rob_6: reg mb2_vdefor i1.y0607##i1.biome##i1.soy_suit i1.post_2007##i1.biome##i1.soy_suit ///
	i1.post_2003##i1.legal_amazon##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)

//
eststo pol_rob_7: reg mb2_vdefor i1.y0607##i1.biome##i.soy_suit i1.post_2007##i1.biome##i.soy_suit ///
	i.post_2007##i1.prodes_mon##i1.soy_suit i1.post_2003##i1.legal_amazon##i1.soy_suit ///
	i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa, vce(cluster municcode)
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
// 	title(Linear regression results\label{tab1}) ///
// 	indicate(Year by state effects = *.year) ///
// 	addnotes("Includes additional covariates and interaction terms as described in SI. Standard errors clustered by municipalities.") ///
	

	
/////// Table S3 - Robustness checks
/// Primary specification	
eststo primary: reg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State by year"
estadd loc cov "Yes"

	
/// No extra covariates
eststo nocov: reg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_state ///
	if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State by year"
estadd loc cov "No"
	
/// Municipality-year fixed effects
eststo muniyear: reg mb2_vdefor i.biome##i.soy_suit##i.post_2005 c.year##i.f_munic ///
	temp trmm roaddist urbandist i.pa i.set if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "Municipality by year"
estadd loc cov "Yes"
	
/// Property-level fixed effects
eststo propfe: xtreg mb2_vdefor i.biome##i.soy_suit##i.post_2005 i.year##i.f_prop_state ///
	temp trmm roaddist urbandist i.pa i.set if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1 & !missing(propid), ///
	fe i(f_propid) vce(cluster f_prop_munic)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State by year"
estadd loc cov "Yes"
	
/// Logistic regression - placeholder for BIFE bias-corrected regression with proper fixed effects from R
eststo logit: logit mb2_vdefor i.biome##i.soy_suit##i.post_2005 ///
	if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode) or
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State and year"
estadd loc cov "Yes"

esttab primary nocov muniyear propfe logit ///
	using `dir'tables\ts4_robustness.tex, replace label se nodepvars fragment ///
	keep(1.biome#1.soy_suit#1.post_2005) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(cov fe n_points N_clust, labels("Includes covariates" "Administrative FE" "N. points" "N. municipalities")) ///
	mtitles("Primary" "No covariates" "Municipality fixed effects" "Property fixed effects" "Logistic") ///
	eform(0 0 0 0 1) eqlabels("" "" "" "" "")

// 	indicate(Covariates = temp) ///
// 	addnotes("Models estimated using forested points located within the Legal Amazon. Standard errors clustered by municipality. Logistic regression (Column 5) depicts odds ratio.") ///

/// Table S4: Robustness to alternate suitability definitions
gen soy_suit2 = soy_suit
eststo prim: reg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez ">Medium"


drop soy_suit2
gen soy_suit2 = suit>0
eststo soares_cat: reg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez "Not included"


drop soy_suit2
gen soy_suit2 = gaezsuit > 40
eststo gaez_cat: reg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares "Not included"
estadd loc gaez ">Medium"

drop soy_suit2
gen soy_suit2 = (suit>0) & (gaezsuit>25)
eststo prim_mod: reg mb2_vdefor i.biome##i.post_2005##i.soy_suit2 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez ">Moderate"

gen soy_suit3 = gaezsuit
eststo gaez_cont: reg mb2_vdefor i.biome##c.soy_suit3##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares "Not included"
estadd loc gaez "Continuous"

drop soy_suit3
gen soy_suit3 = (suit>0) * gaezsuit
eststo both_cont: reg mb2_vdefor i.biome##c.soy_suit3##i.post_2005 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc soares ">Suitable"
estadd loc gaez "Continuous"

label var soy_suit2 "Suitable for soy"
label values soy_suit2 soy_labels
label var soy_suit3 "Continuous suitability"

esttab prim soares_cat gaez_cat prim_mod gaez_cont both_cont ///
	using `dir'tables\ts5_soy_robustness.tex, replace label wrap nomtitles fragment ///
	keep(1.biome#1.post_2005#1.soy_suit2 1.biome#1.post_2005#c.soy_suit3) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) modelwidth(8) ///
	stats(soares gaez n_points N_clust, labels("Soares-Filho et al data treatment" "GAEZ data treatment" "N. points" "N. municipalities"))
// 	addnotes("Models estimated using forested points located within the Legal Amazon. Standard errors clustered by municipality.") ///
// 	mtitles("Primary" "No covariates" "Municipality fixed effects" "Property fixed effects" "Logistic") ///
// 	eform(0 0 0 0 1)

drop soy_suit2 soy_suit3

////// Figure 2C: Time effects
eststo time: reg mb2_vdefor b2005.year##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa  if legal_amazon==1 & dist_amb>-300 & dist_aml<300, vce(cluster municcode)	
esttab time ///
	using `dir'figures\f2_time_plot.csv, replace plain ///
	b(a3) ci(a3)



