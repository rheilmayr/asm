clear all
clear matrix
clear mata
set matsize 11000
set maxvar 10000
eststo clear
set more off

loc dir D:\cloud\Dropbox\collaborations\glue-sb\soyM\analysis\public\
insheet using `dir'long.csv

drop if dist_amb<-300
drop if dist_aml>300
drop if year<=2001

gen post_2004 = year>2004
gen post_2006 = year>2006
gen post_2008 = year>2008
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
label var post_2006 "Post-2006"
label var post_2004 "Post-2004"
label define biome_labels 0 "Cerrado biome" 1 "Amazon biome"
label values biome biome_labels
label define l_labels 0 "Not legal Amazon" 1 "Legal Amazon"
label values legal_amazon l_labels
label define p_labels 0 "Pre-2006" 1 "Post-2006"
label values post_2006 p_labels
label define p4_labels 0 "Pre-2004" 1 "Post-2004"
label values post_2004 p4_labels
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

/// Table 1: Triple diffs
eststo clear
eststo dd_ss_inbio: reg mb2_vdefor i.soy_suit##i.post_2006 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & biome==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Within Amazon biome"
	
eststo dd_ss_outbio: reg mb2_vdefor i.soy_suit##i.post_2006 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & biome==0, ///
	vce(cluster municcode)		
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"
	
eststo dd_bio_inss: reg mb2_vdefor i.biome##i.post_2006 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & soy_suit==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Within soy-suitable"
	
eststo dd_bio_outss: reg mb2_vdefor i.biome##i.post_2006 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1 & soy_suit==0, ///
	vce(cluster municcode)		
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside soy-suitable"

eststo ddd: reg mb2_vdefor i.biome##i.soy_suit##i.post_2006 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "All points"
	
// Export table
esttab using `dir'tables\t1_ddd.tex, replace label nodepvars ///
	keep(1.soy_suit#1.post_2006 1.biome#1.post_2006 1.biome#1.soy_suit#1.post_2006) ///
	booktabs width(0.8\hsize) alignment(c) ///
	title(Linear regression results\label{tab1}) ///
	indicate(Year by state effects = *.year) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(sample n_points N_clust, labels("Sample" "N. points" "N. municipalities")) ///
	addnotes("Includes additional covariates and interaction terms as described in SI. Standard errors clustered by municipality") ///
	mtitles("Amazon DD" "Cerrado DD" "Soy-suitable DD" "Non-soy suitable DD" "Triple difference")
eststo clear

//// Table 2: Policy interactions
gen on_car = !missing(propid)
gen on_car2 = !missing(car_year)
replace on_car2 = . if (state!="MT" & state!="PA")
gen asm_now = post_2006 * biome

scalar drop _all
program drop _all 
program define prog1, rclass
	scalar p = 2*ttail(r(df),abs(r(estimate)/r(se)))
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
	return local coef `=scalar(coef)'
end	

reg mb2_vdefor i.asm_now i.biome i.year ///
	temp trmm urbandist roaddist i.pa i.set if soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo r0
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_t `r(t)'

reg mb2_vdefor i.asm_now##i.car_now i.biome##i.on_car2 i.year ///
	temp trmm urbandist roaddist i.pa i.set if soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo r1
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_t `r(t)'
lincom 1.car_now
prog1
estadd local c `r(coef)'
estadd local c_t `r(t)'
lincom 1.asm_now + 1.car_now + 1.asm_now#1.car_now
prog1
estadd local ac `r(coef)'
estadd local ac_t `r(t)'

reg mb2_vdefor i.asm_now##i.gts_now i.biome##i.gts_ever i.year##i.f_state ///
	temp trmm urbandist roaddist i.pa i.set if year>2001 & soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo r2
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_t `r(t)'
lincom 1.asm_now + 1.gts_now + 1.asm_now#1.gts_now
prog1
estadd local ag `r(coef)'
estadd local ag_t `r(t)'

reg mb2_vdefor i.asm_now##i.gts_now##i.car_now i.biome##i.gts_ever##i.on_car2 i.year##i.f_state ///
	temp trmm urbandist roaddist i.pa i.set if year>2001 & soy_suit==1 & legal_amazon==1 & (state=="MT" | state=="PA"), vce(cluster municcode)
count if (year==2002 & e(sample)==1)
	estadd scalar n_points = r(N)
eststo r3
lincom 1.asm_now
prog1
estadd local a `r(coef)'
estadd local a_t `r(t)'
lincom 1.car_now
prog1
estadd local c `r(coef)'
estadd local c_t `r(t)'
lincom 1.asm_now + 1.gts_now + 1.asm_now#1.gts_now
prog1
estadd local ag `r(coef)'
estadd local ag_t `r(t)'
lincom 1.asm_now + 1.car_now + 1.asm_now#1.car_now
prog1
estadd local ac `r(coef)'
estadd local ac_t `r(t)'
lincom 1.asm_now + 1.gts_now + 1.car_now + 1.asm_now#1.car_now + 1.gts_now#1.car_now + 1.asm_now#1.gts_now#1.car_now + 1.asm_now#1.gts_now
prog1
estadd local acg `r(coef)'
estadd local acg_t `r(t)'

// Export table
esttab using `dir'tables\t2_complementary_policies.tex, replace label nodepvars nomtitles ///
	keep(1.asm_now) ///
	booktabs width(0.8\hsize) alignment(c) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	addnotes("Models estimated using forested, soy-suitable points in the states of Mato Grosso and Para. Models include additional covariates and interaction terms as described in SI. Standard errors were clustered by municipalities." ///
	"Terms presented in this table represent linear combinations of coefficiencts on individual and interaction terms to capture aggregate effect of multiple policies.") ///
	stats(a a_t c c_t ac ac_t ag ag_t acg acg_t n_points N_clust, ///
	layout(@ (@) @ (@) @ (@) @ (@) @ (@) @ @) ///
	labels("ASM only" " " "CAR only" " " "ASM and CAR" " " "ASM and GTS" " " "ASM, CAR and GTS" " " ///
	"N. points" "N. municipalities"))

///////// Table 3: Leakage	
eststo clear
// ILUC
reg mb2_vdefor i.post_2006##i.biome i.year##i.f_state ///
	if legal_amazon==1 & soy_suit==0, ///
	vce(cluster municcode)
eststo r1
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside soy-suitable"

// Outside of gts monitored areas
reg mb2_vdefor i.post_2006##i.biome i.year##i.f_state ///
	if legal_amazon==1 & soy_suit==1 & gts_ever==0, ///
	vce(cluster municcode)	
eststo r3
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Soy-suitable, outside GTS"
	
// traditional distance leakage analysis
reg mb2_vdefor c.proximity##i.post_2006 i.year##i.f_state ///
	if legal_amazon==1 & biome==0 & soy_suit==1, ///
	vce(cluster municcode)
eststo r4
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"
	
reg mb2_vdefor i.post_2006##i.close i.year##i.f_state ///
	if biome==0 & legal_amazon==1 & soy_suit==1, ///
	vce(cluster municcode)	
eststo r6	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample "Outside Amazon biome"

esttab using `dir'tables\t3_leak.tex, replace label nodepvars nomtitles ///
	keep(1.post_2006#1.biome 1.post_2006#c.proximity 1.post_2006#1.close) ///
	booktabs width(0.8\hsize) alignment(c) ///
	title(Linear regression results\label{tab1}) ///
	indicate(Year by state effects = *.year) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(sample n_points N_clust, labels("Sample" "N. points" "N. municipalities")) ///
	addnotes("Models estimated using forested points located within the Legal Amazon. Models include additional covariates and interaction terms as described in SI. Standard errors clustered by municipality.")
eststo clear	
		
	
///////// Table S1: Summary stats
eststo clear
reg mb2_vdefor i.biome##i.soy_suit##i.post_2006 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if legal_amazon==1, ///
	vce(cluster municcode)
estpost summarize mb2_vdefor biome legal_amazon soy_suit temp trmm roaddist urbandist set pa
esttab using `dir'tables\ts1_summary.tex, cells("mean(fmt(a3)) sd(fmt(a3)) min(fmt(a3)) max(fmt(a3))") nomtitle nonumber ///
	varlabels(mb2_vdefor "Deforested (binary)" biome "In Amazon biome (binary)" legal_amazon "In Legal Amazon (binary)" ///
	soy_suit "Suitable for soy (binary)" temp "Mean temperature (degrees C)" trmm "Annual precipitation (cm/y)" ///
	roaddist "Distance to nearest road (km)" urbandist "Distance to nearest town (km)" set "In settlement (binary)" ///
	pa "In protected area (binary)") labcol2("Mapbiomas" "MMA, 2018" "MMA, 2018" ///
	"Soares-Filho et al., 2014" "MODIS" "TRMM" "IBGE, 2010" "DNIT, 2018" ///
	"INCRA, 2018" "Funai, 2018 and MMA, 2018", title(Source)) stats() replace
	
	
////// Table S2: Legal amazon and beyond
// Within legal amazon, homogeneous effect
eststo lam_homog: reg mb2_vdefor i1.post_2006##i1.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if legal_amazon==1, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample = "Within Legal Amazon"
	
// Within legal amazon, heterogeneous effect
eststo lam_hetero: reg mb2_vdefor i1.post_2006##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa if legal_amazon==1, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample = "Within Legal Amazon"

// All points, homogeneous effect
eststo all_homog: reg mb2_vdefor i1.post_2006##i1.biome i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample = "All points"
	
// All points, heterogeneous effect
eststo all_hetero: reg mb2_vdefor i1.post_2006##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample = "All points"

// Legal amazon interaction, homogeneous effect	
eststo lamint_homog: reg mb2_vdefor i1.post_2006#i1.biome i1.biome ///
	i1.post_2004#i1.legal_amazon i1.legal_amazon i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample = "All points"
	
// Legal amazon interaction, heterogeneous effect
eststo lamint_hetero: reg mb2_vdefor i1.post_2006#i1.biome i1.post_2006#i1.biome#i1.soy_suit ///
	i1.post_2004#i1.legal_amazon i1.post_2004#i1.legal_amazon#i1.soy_suit ///
	i1.biome#i1.soy_suit i1.legal_amazon#i1.soy_suit i1.post_2006#i1.soy_suit i1.post_2004#i1.soy_suit ///
	i1.biome i1.soy_suit i1.legal_amazon i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa, vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc sample = "All points"

// Export table
esttab using `dir'tables\ts2_ppcdam_interact.tex, replace label nodepvars nomtitles ///
	keep(1.post_2006#1.biome 1.post_2006#1.biome#1.soy_suit ///
		 1.post_2004#1.legal_amazon 1.post_2004#1.legal_amazon#1.soy_suit) ///
	order(1.post_2006#1.biome 1.post_2006#1.biome#1.soy_suit ///
		1.post_2004#1.legal_amazon 1.post_2004#1.legal_amazon#1.soy_suit) ///
	booktabs width(0.8\hsize) alignment(c) ///
	title(Linear regression results\label{tab1}) ///
	indicate(Year by state effects = *.year) ///
	stats(n_points N_clust, labels("N. points" "N. municipalities")) ///
	addnotes("Includes additional covariates and interaction terms as described in SI. Standard errors clustered by municipalities.") ///
	mgroups("Only Legal Amazon" "All points" "All points with interactions", pattern(1 0 1 0 1 0) prefix(\multicolumn{@span}{c}{) suffix(}) span erepeat(\cmidrule(lr){@span})) ///
	star(* 0.1 ** 0.05 *** 0.01)
	
eststo clear

	
/////// Table S3 - Robustness checks
eststo clear
/// Primary specification	
eststo primary: reg mb2_vdefor i.biome##i.soy_suit##i.post_2006 i.year##i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State by year"
estadd loc cov "Yes"

	
/// No extra covariates
eststo nocov: reg mb2_vdefor i.biome##i.soy_suit##i.post_2006 i.year##i.f_state ///
	if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State by year"
estadd loc cov "No"
	
/// Municipality-year fixed effects
eststo muniyear: reg mb2_vdefor i.biome##i.soy_suit##i.post_2006 c.year##i.f_munic ///
	temp trmm roaddist urbandist i.pa i.set if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode)	
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "Municipality by year"
estadd loc cov "Yes"
	
/// Property-level fixed effects
eststo propfe: xtreg mb2_vdefor i.biome##i.soy_suit##i.post_2006 i.year##i.f_prop_state ///
	temp trmm roaddist urbandist i.pa i.set if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1 & !missing(propid), ///
	fe i(f_propid) vce(cluster f_prop_munic)
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State by year"
estadd loc cov "Yes"
	
/// Logistic regression
eststo logit: logit mb2_vdefor i.biome##i.soy_suit##i.post_2006 i.year i.f_state ///
	temp trmm roaddist urbandist i.pa i.set if year>2001 & dist_amb>-300 & dist_aml<300 & legal_amazon==1, ///
	vce(cluster municcode) or
count if (year==2002 & e(sample)==1)
estadd scalar n_points = r(N)
estadd loc fe "State and year"
estadd loc cov "Yes"

esttab using `dir'tables\ts4_robustness.tex, replace label ///
	keep(1.biome#1.soy_suit#1.post_2006) ///
	booktabs width(0.8\hsize) alignment(c) ///
	indicate(Covariates = temp) ///
	star(* 0.1 ** 0.05 *** 0.01) ///
	stats(fe n_points N_clust, labels("Administrative FE" "N. points" "N. municipalities")) ///
	addnotes("Models estimated using forested points located within the Legal Amazon. Standard errors clustered by municipality. Logistic regression (Column 5) depicts odds ratio.") ///
	mtitles("Primary" "No covariates" "Municipality fixed effects" "Property fixed effects" "Logistic") ///
	eform(0 0 0 0 1)

	
////// Figure 2C: Time effects
eststo clear
eststo time: reg mb2_vdefor b2006.year##i1.biome##i1.soy_suit i.year##i.f_state ///
	temp trmm roaddist urbandist i.set i.pa  if legal_amazon==1 & dist_amb>-300 & dist_aml<300, vce(cluster municcode)	
esttab using `dir'figures\f2_time_plot.csv, replace plain ///
	b(a3) ci(a3)
