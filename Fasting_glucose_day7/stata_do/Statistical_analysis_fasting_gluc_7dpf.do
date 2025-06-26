set more off
clear
cd "/Users/marde358/Desktop/Projects/Novo_screen"

foreach x in "SIRT1" "ADRA2A" "NOTCH2" "SLC2A2" "POLDIP2" ///
"ACRV1C" "TBC1D4" "TP53INP1" "DUSP26" ///
"IGF1" "SLC2A10" "ATP2A3" "BCL11A" "ZNF598" "LMNA" ///
"GCG" "IRX3" "GLP1R" "PLIN1" ///
"SLC5A2" "LPL" "PDX1" "GCK" ///
{

** Examine the effect of mutations in each gene on glucose adjusted for batch (order in which the plate was read) and injection set (if it was injected in two different clutches) 

use "DATA/glucose_day7/DATA_glucose_bygene/`x'/`x'_for_analysis.dta", clear

gen x="`x'"
*egen y=ends(x) , punct(_) tail trim
*levelsof y, local(y) clean
gen z=lower("`x'")
levelsof z, local(z) clean

foreach outcome in gluc intgluc {

preserve
regress `outcome' i.`z'_dich batch i.injset position
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/glucose_day7/`x'/`z'_dich_batch/`outcome'_est", replace)
use "RESULTS/glucose_day7/`x'/`z'_dich_batch/`outcome'_est", clear
gen N=e(N)
save "RESULTS/glucose_day7/`x'/`z'_dich_batch/`outcome'_est", replace
restore
}

* now append
clear

save "RESULTS/glucose_day7/`x'/`z'_dich_batch/`z'_dich_batch.dta", replace emptyok 

foreach file in gluc_est intgluc_est  {

use "RESULTS/glucose_day7/`x'/`z'_dich_batch/`file'", clear
gen model="`file'"

append using "RESULTS/glucose_day7/`x'/`z'_dich_batch/`z'_dich_batch.dta"
save "RESULTS/glucose_day7/`x'/`z'_dich_batch/`z'_dich_batch.dta", replace
}
use "RESULTS/glucose_day7/`x'/`z'_dich_batch/`z'_dich_batch.dta", clear
drop if strpos(parm,"b.")
drop if strpos(parm,"o.")
drop dof
egen id=concat(parm) if parm=="_cons"
replace id=parm if id==""
sort model, stable
by model: gen order=_n

gen trait = substr(model,1,length(model)-4)
rename t z
order trait model parm estimate stderr min95 max95 t p N
rename (estimate stderr min95 max95 N p z) (beta_batch SE_batch LCI_batch UCI_batch N_batch P_batch z_batch)
sort trait model order
drop id order
save "RESULTS/glucose_day7/`x'/fixed_effects_f0_`z'_dich_batch.dta", replace
export excel using "RESULTS/glucose_day7/`x'/fixed_effects_f0_`z'_dich_batch.xls", replace firstrow(variables)



** Merge results files with summary statistics
use "RESULTS/glucose_day7/`x'/fixed_effects_f0_`z'_dich_batch.dta", clear
*erge 1:1 trait model parm using "RESULTS/liver/ins_lfabp10a/`x'/fixed_effects_f0_`y'_dich_size.dta"
*drop _merge

gen order = .
replace order = 1 if strpos(parm,"1.`z'_dich")
replace order = 2 if strpos(parm,"batch")
replace order = 3 if strpos(parm,"_cons")

gen trait_order = .

replace trait_order = 1 if strpos(trait,"gluc")
replace trait_order = 2 if strpos(trait,"intgluc")

replace parm = "intercept" if parm == "_cons"

* make some changes to prevent work on results table

foreach var in beta_batch SE_batch LCI_batch UCI_batch z_batch P_batch N_batch   {
replace `var´' = . if parm == "intercept"
}

foreach var in trait model parm {
replace `var´' = "" if parm == "intercept"
}

gen l=""

sort trait_order model order parm
drop order trait_order

format beta_batch SE_batch LCI_batch UCI_batch z_batch    %6.3f
format P_batch  %12.2e
drop if model == "random" & parm == ""
drop model
drop if beta_batch == .
order trait parm beta_batch SE_batch LCI_batch UCI_batch z_batch P_batch N_batch l     
replace parm = "`z'_dich" if parm == "1.`z'_dich"
*drop if strpos(parm,".batch")
*drop if strpos(parm,".tank")
*drop eq
drop l
* save sumstats 
outsheet using "RESULTS/glucose_day7/`x'/`z'_dich_results.txt", noquote replace
}




