set more off
clear
cd "C:\Users\endmu284\Documents\Size\Data\S026_LPL\"


*create a label for conditions in order to plot
use "VAST_cleargroup.dta", clear
*drop conditions
gen conditions = "Controls" if condition == 0
replace conditions = "lpl_crispants" if condition == 1

save "VAST_cleargroup.dta", replace


** Examine the effect of mutants vs kita  on outcomes

use "VAST_cleargroup.dta", clear

foreach outcome in intlength ///
intangle 
 {


preserve
regress `outcome' i.condition tod 
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS\S026_LPL`outcome'_est", replace)
use "RESULTS\S026_LPL`outcome'_est", clear
gen N=e(N)
save "RESULTS\S026_LPL`outcome'_est.dta", replace
restore
}


* now append
clear
save "RESULTS\S026_LPL.dta", replace emptyok 

foreach file in intlength_est ///
intangle_est 
{

use "RESULTS\S026_LPL`file'", clear
gen model="`file'"
append using "RESULTS\S026_LPL.dta"
save "RESULTS\S026_LPL.dta", replace
}
use "RESULTS\S026_LPL.dta", clear
drop if strpos(parm,"b.")
drop if strpos(parm,"o.")
drop dof
egen id=concat(parm) if parm=="_cons"
replace id=parm if id==""
sort model, stable
by model: gen order=_n


gen  trait = substr(model,4,length(model)-4)
order trait model parm estimate stderr min95 max95 t p N
rename (estimate stderr min95 max95 p) (beta SE LCI UCI P)
sort trait model order
drop id order


gen trait_order = .

replace trait_order = 1 if strpos(trait,"length")
replace trait_order = 2 if strpos(trait,"angle")


replace parm = "lpl_crispants" if parm == "1.condition"
replace parm = "tod" if parm == "tod"


* make some changes to prevent work on results table
replace parm = "intercept" if parm == "_cons"

foreach var in beta SE LCI UCI t P N {
replace `var´' = . if parm == "intercept"
}

foreach var in trait model parm {
replace `var´' = "" if parm == "intercept"
}

sort trait_order model  parm
drop trait_order

format beta SE LCI UCI t %6.3f
format P %12.2e
*drop if model == "random" & parm == ""
*drop model
order trait parm beta SE LCI UCI t P N
drop if beta == .
*replace trait = substr(trait,4,.) if strpos(trait,"int")


* save sum stats for Novo and for figures
save "RESULTS\S026_LPL_results.dta", replace
outsheet using "RESULTS\S026_LPL_results.txt", noquote replace
outsheet using "C:\Users\endmu284\Documents\Size\Data\RESULTS\S026_LPL_results.txt", noquote replace


