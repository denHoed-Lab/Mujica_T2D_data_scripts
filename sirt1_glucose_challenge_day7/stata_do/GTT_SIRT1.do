cd "PASTE MAIN PATH HERE"


* run the analysis for Glucose Tolernance Test for SIRT1 gene and generate figure 3. 
clear
save "something.dta", emptyok replace

foreach x in  "2_SIRT1"  {
use "something.dta", clear
set obs 1
gen x="`x'"
egen w=ends(x) , punct(_) tail trim
levelsof w, local(w) clean
gen y=lower("`w'")
levelsof y, local(y) clean

use "`x'/`y'_forGTT.dta", clear

*Generate position variable for the time/position in the 384 well plate
gsort  row384 column384
gene batch1 = _n
gsort  row384 -column384
gene batch2 = _n
tab row384
replace batch1= 0 if inlist(row384, "B", "D", "F", "H", "J", "L", "N", "P")
replace batch2= 0 if !inlist(row384, "B", "D", "F", "H", "J", "L", "N", "P")
gen position = batch1+batch2
gsort position
drop batch1 batch2


local ytitle_m = "Glucose (µM/sample)"
local ytitle_c = "Contrast of linear prediction"

drop if time==0
replace time =time -60

regress gluc time##i.`y'_dich i.assay_plate position
margins time#`y'_dich, saving("`y'_margins.dta", replace)

preserve
use "`y'_margins.dta", clear
rename (_margin _ci_lb _ci_ub _m1 _m2) (margin LCI UCI time `y'_dich)
drop _deriv _term _predict _at _atopt _se _statistic _pvalue

replace time = time -2 if `y'_dich == 0
replace time = time +2 if `y'_dich == 1

* marginsplot
twoway scatter margin time if `y'_dich==0, connect(line) mcolor(gray) ms(Oh) lpattern(shortdash) lcolor(gray) xlabel(0 60 90 120 180) || ///
rcap LCI UCI time if `y'_dich==0, lcolor(gray) || ///
scatter margin time if `y'_dich==1, connect(line) mcolor(orange) lpattern(shortdash) lcolor(orange) || ///
rcap LCI UCI time if `y'_dich==1, lcolor(orange) ///
name(mp_`y', replace) graphregion(color(white)) ///
legend(off) ylabel(, nogrid) xlabel(, nogrid) /// 
ytitle("`ytitle_m'") xtitle(" ") yscale(range(10 17)) ylabel(5 (5) 17)
restore



* contrast crispants and controls
margins r.`y'_dich@time, contrast(nowald) saving("`y'_contrast.dta", replace)

preserve
use "`y'_contrast.dta", clear
rename (_margin _ci_lb _ci_ub _m2) (margin LCI UCI time)
drop _deriv _term _predict _at _atopt _se _statistic _pvalue _m1

* contrast plot
twoway scatter margin time, connect(line) mcolor(orange)  lpattern(shortdash) lcolor(orange) xlabel(0 60 90 120 180) || ///
rcap LCI UCI time, lcolor(orange) ///
name(contrast_`y', replace) graphregion(color(white)) ///
legend(off) ylabel(, nogrid) xlabel(, nogrid) yline(0, lcolor(gray)) /// 
ytitle("`ytitle_c'") xtitle("time") yscale(range(-5 10)) ylabel(-5 (5) 10)
restore

* Graph combine
graph combine mp_`y' contrast_`y', cols(1) xsize(14) ysize(14) graphregion(color(white) lwidth(large)) title({it:`y'}) name("`y'_uM", replace)

graph export "/Users/endmu284/Desktop/Projects/Novo_screen/FIGURES_GTT/results_GTT_`y'_uM.png", width(2000) height(2000) replace
graph export "/Users/endmu284/Desktop/Projects/Novo_screen/FIGURES_GTT/results_GTT_`y'_uM.svg", width(2000) height(2000) replace
}



* run the analysis
clear
save "something.dta", emptyok replace

foreach x in "2_SIRT1"  {
use "something.dta", clear
set obs 1
gen x="`x'"
egen w=ends(x) , punct(_) tail trim
levelsof w, local(w) clean
gen y=lower("`w'")
levelsof y, local(y) clean

use "`x'/`y'_forGTT.dta", clear

*Generate position variable for the time/position in the 384 well plate
gsort  row384 column384
gene batch1 = _n
gsort  row384 -column384
gene batch2 = _n
tab row384
replace batch1= 0 if inlist(row384, "B", "D", "F", "H", "J", "L", "N", "P")
replace batch2= 0 if !inlist(row384, "B", "D", "F", "H", "J", "L", "N", "P")
gen position = batch1+batch2
gsort position
drop batch1 batch2


local ytitle_m = "Glucose (µM/sample)"
local ytitle_c = "Contrast of linear prediction"

regress gluc time#i.`y'_dich i.assay_plate position
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e es_1) saving ("`x'_est_180", replace)
use "`x'_est_180", clear
gen N=e(N)
save "`x'_est", replace
export excel using "`x'_est_180.xls", replace firstrow(variables)


}
