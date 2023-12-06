*********************
* Determinants and heterogeneity of switching costs in IT outsourcing: estimates from firm-level data *
* Christian Peukert *
* European Journal of Information Systems, 2018 *
* https://doi.org/10.1080/0960085X.2018.1529374 *
*********************

* This code reproduces the results reported in Table 2, 3 and 5

set more off
version 10
global texpath="results"

use "data.dta",clear

qui probit $v px $z $ind $other
reg $c px cx $z $other if e(sample), nocons
estimates store c

version 10
probit $v px $z $ind $other if e(sample)
utest lvex lvex2, fieller
utest swex swex2, fieller
version 11

probit $v px $z $ind $other if e(sample)
global k=e(k)
estimates store v

gen smpl=e(sample)
suest c
estimates store crobust
suest v
estimates store vrobust
suest c v
matrix V=e(V)
matrix b=e(b)
erepost b=b V=V
estimates store sandwich

local alpha="[c_mean]px+[c_mean]cx"
local sigma="(`alpha')/[v_${v}]px"

	qui {
		global nlcom "(alpha: `alpha') (sigma: `sigma')"
		foreach ind in $z $ind {
			global nlcom "${nlcom} (`ind': `sigma'*[v_${v}]`ind')"
		}
	}
global nlcom "$nlcom (_cons: `sigma'*[v_${v}]_cons)"
nlcom ${nlcom},post
estimates store nlcom_sandwich


/// Tables 2, 3

//--------------------------------------------
// Standard error clustering



qui { 
	suest c v,vce(cluster join_number)
	matrix V_1=e(V)
	suest c v,vce(cluster vidt)
	matrix V_2=e(V)
	suest c v,vce(cluster join_number_vidt)
	matrix V_12=e(V)
	matrix b=e(b)
	matrix V=V_1+V_2-V_12
	erepost b=b V=V
	estimates store cluster
}

local alpha="[c_mean]px+[c_mean]cx"
local sigma="(`alpha')/[v_${v}]px"


nlcom ${nlcom},post
estimates store nlcom_cluster
	matrix B=e(b)

//--------------------------------------------
// Standardized coefficients in order to be able to compare coefficients

qui { 
	suest c v,vce(cluster join_number)
	matrix V_1=e(V)
	suest c v,vce(cluster vidt)
	matrix V_2=e(V)
	suest c v,vce(cluster join_number_vidt)
	matrix V_12=e(V)
	matrix b=e(b)
	matrix V=V_1+V_2-V_12
	erepost b=b V=V
	estimates store cluster
}

local alpha="[c_mean]px+[c_mean]cx"
local sigma="(`alpha')/[v_${v}]px"

qui{
foreach var in $z $ind {
	qui su `var' if smpl
	local `var'_sd=r(sd)
	local `var'_dummy=r(min)==0 & r(max)==1
}
}

global nlcom_std "(_cons: (`sigma')*[v_${v}]_cons)"
global nlcom_std "${nlcom_std} (llogbarg: 2*`llogbarg_sd'*`sigma'*[v_${v}]llogbarg)"
global nlcom_std "${nlcom_std} (hhi_cont: 2*`hhi_cont_sd'*`sigma'*[v_${v}]hhi_cont)"
global nlcom_std "${nlcom_std} (dvrepd: `dvrepd_sd'*`sigma'*[v_${v}]dvrepd)"
global nlcom_std "${nlcom_std} (dnostd: 2*`dnostd_sd'*`sigma'*[v_${v}]dnostd)"
global nlcom_std "${nlcom_std} (logemp: 2*`logemp_sd'*`sigma'*[v_${v}]logemp)"
global nlcom_std "${nlcom_std} (lhhi: 2*`lhhi_sd'*`sigma'*[v_${v}]lhhi)"
global nlcom_std "${nlcom_std} (lspec: 2*`lspec_sd'*`sigma'*[v_${v}]lspec)"
global nlcom_std "${nlcom_std} (lvex_s: 2*`lvex_sd'*`sigma'*[v_${v}]lvex+2*`lvex_sd'^2*`sigma'*[v_${v}]lvex2)"
global nlcom_std "${nlcom_std} (logtravelemp: 2*`logtravelemp_sd'*`sigma'*[v_${v}]logtravelemp)"
global nlcom_std "${nlcom_std} (swex_s: 2*`swex_sd'*`sigma'*[v_${v}]swex+2*`swex_sd'^2*`sigma'*[v_${v}]swex2)"
		

nlcom ${nlcom_std}, post
estimates store nlcom_cluster_std


esttab c v sandwich cluster using "table_2.tex", replace wide noomitted order(px cx $order) mtitle("Sandwich" "Sandwich" "Sandwich" "Cluster") indicate("Year fixed effect=y*" "Vendor fixed effect t=vmatch*" "Vendor fixed effect t-1=lvmatch*") se star(* 0.10 ** 0.05 *** 0.01) scalars("r2_a adj. R2" "r2_p pseudo R2") obslast label compress booktabs equations(1:.:1:1,.:1:3:3) eqlabels(c v)
esttab nlcom_cluster nlcom_cluster_std using "table_3.tex", replace wide order(alpha sigma $order) mtitle("Cluster" "Cluster Std") se star(* 0.10 ** 0.05 *** 0.01) obslast label compress booktabs

/// Table 5


gen s=0

local i=0
foreach var in 0 0 0 0 0 0 0 $s 0{
	local i=`i'+1
	di B[1,`i'] "`var'"
	replace s=s-B[1,`i']*`var' if smpl
}

gen srel=s/mean_acct_290d
gen srelt=s/mean_acct_290d/vex
gen srelpx=s/px


label variable s "Absolute"
label variable srel "Relative"
label variable srelt "Relative per Contract Year"
label variable srelpx "Relative $\p_\tau$"
label define vswm 0 "Non-Switchers" 1 "Switchers"
label values vswm vswm

estpost tabstat s srel srelt, stats(mean sd min max p50 skewness kurtosis iqr) by(vswm) columns(s)
esttab using "table_5.tex", replace ////
cells("mean(fmt(%9.4g) label(Mean)) sd(fmt(%9.4g) label(SD)) min(fmt(%9.4g) label(Min)) max(fmt(%9.4g) label(Max)) p50(fmt(%9.4g) label(Median)) iqr(fmt(%9.4g) label(IQR))") ////
nomtitle nonumber booktabs label noobs
