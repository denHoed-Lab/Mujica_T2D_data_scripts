set more off
clear
cd "Paste_Main_Path_Here"

* data for these genes are analysed separately because for a subset of larvae targeted larvae and controls were raised in separate tanks
* as a result, the analysis should not be adjusted for tank for these genes.
 
foreach x in "23_TP53INP1" "33_POLDIP2" "44_ANXA6" "127_PDX1" "150_GCK" {
	
** Examine the effect of mutations in each gene adjusted for tank, batch (in residuals) and tod
use "DATA/ins_lfabp10a_bygene/`x'/`x'_for_analysis.dta", clear

gen x="`x'"
egen y=ends(x) , punct(_) tail trim
levelsof y, local(y) clean
gen z=lower("`y'")
levelsof z, local(z) clean

foreach outcome in intlength intprotein ///
intangle inteye_area ///
intLiverArea intLiverLipN intLiverLipArea ///
intNrCells intVolume intCellVolume_Av intCellInt_Tot intCellInt_Av intIV intISA intIED intIPAL1 intIPAL2 intIPAL3 {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace
restore
}

foreach outcome in intldl inttg inttc intgluc {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age pos_bio i.batch_bio
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace
restore
}
 
foreach outcome in intdorsal_area intlateral_area {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age length
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace
restore
}

foreach outcome in intLiverLipInt_Tot {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age LiverLipIntBG
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace
restore
}

foreach outcome in intpig_dorsal_tail {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age dorsal_tail
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace
restore
}

foreach outcome in vld_area {

preserve
capture nbreg `outcome' i.`z'_dich tod i.batch i.age, irr
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`outcome'_est", replace
restore
}

* now append
clear
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`z'_dich_batch.dta", replace emptyok 

foreach file in intlength_est intdorsal_area_est intlateral_area_est intprotein_est ///
intangle_est inteye_area_est intpig_dorsal_tail_est ///
intldl_est inttg_est inttc_est intgluc_est ///
intLiverArea_est intLiverLipN_est intLiverLipArea_est intLiverLipInt_Tot_est ///
intNrCells_est intVolume_est intCellVolume_Av_est intCellInt_Tot_est intCellInt_Av_est intIV_est intISA_est intIED_est intIPAL1_est intIPAL2_est intIPAL3_est ///
vld_area_est {

use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`file'", clear
gen model="`file'"
append using "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`z'_dich_batch.dta"
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`z'_dich_batch.dta", replace
}
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_batch/`z'_dich_batch.dta", clear
drop if strpos(parm,"b.")
drop if strpos(parm,"o.")
drop dof
egen id=concat(parm) if parm=="_cons"
replace id=parm if id==""
sort model, stable
by model: gen order=_n

gen trait = substr(model,1,length(model)-4)
order trait model parm estimate stderr min95 max95 t p N
rename (estimate stderr min95 max95 N p t) (beta_batch SE_batch LCI_batch UCI_batch N_batch P_batch t_batch)
sort trait model order
drop id order
save "RESULTS/ins_lfabp10a/`x'/fixed_effects_f0_`z'_dich_batch.dta", replace
export excel using "RESULTS/ins_lfabp10a/`x'/fixed_effects_f0_`z'_dich_batch.xls", replace firstrow(variables)





******************************





** Examine the effect of mutations additionally adjusted for size
use "DATA/ins_lfabp10a_bygene/`x'/`x'_for_analysis.dta", clear

gen x="`x'"
egen y=ends(x) , punct(_) tail trim
levelsof y, local(y) clean
gen z=lower("`y'")
levelsof z, local(z) clean

foreach outcome in intLiverArea intLiverLipN intLiverLipArea ///
intNrCells intVolume intCellVolume_Av intCellInt_Tot intCellInt_Av intIV intISA intIED intIPAL1 intIPAL2 intIPAL3 {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age length dorsal_area
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace
restore
}

foreach outcome in intldl inttg inttc intgluc {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age length dorsal_area pos_bio i.batch_bio
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace
restore
}
 
foreach outcome in intLiverLipInt_Tot {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age LiverLipIntBG length dorsal_area
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace
restore
}

foreach outcome in intpig_dorsal_tail {

preserve
regress `outcome' i.`z'_dich tod i.batch i.age dorsal_tail length dorsal_area
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace
restore
}

foreach outcome in vld_area {

preserve
capture nbreg `outcome' i.`z'_dich tod i.batch i.age length dorsal_area, irr
parmest, format(estimate %8.6f min95 max95 %8.6f p %8.1e) saving ("RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace)
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", clear
gen N=e(N)
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`outcome'_est", replace
restore
}

* now append
clear
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`z'_dich_size.dta", replace emptyok 

foreach file in intldl_est inttg_est inttc_est intgluc_est ///
intLiverArea_est intLiverLipN_est intLiverLipArea_est intLiverLipInt_Tot_est ///
intNrCells_est intVolume_est intCellVolume_Av_est intCellInt_Tot_est intCellInt_Av_est intIV_est intISA_est intIED_est intIPAL1_est intIPAL2_est intIPAL3_est ///
vld_area_est {

use "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`file'", clear
gen model="`file'"
append using "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`z'_dich_size.dta"
save "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`z'_dich_size.dta", replace
}
use "RESULTS/ins_lfabp10a/`x'/`z'_dich_size/`z'_dich_size.dta", clear
drop if strpos(parm,"b.")
drop if strpos(parm,"o.")
drop dof
egen id=concat(parm) if parm=="_cons"
replace id=parm if id==""
sort model, stable
by model: gen order=_n

gen trait = substr(model,1,length(model)-4)
order trait model parm estimate stderr min95 max95 t p N
rename (estimate stderr min95 max95 N p t) (beta_size SE_size LCI_size UCI_size N_size P_size t_size)
sort trait model order
drop id order
save "RESULTS/ins_lfabp10a/`x'/fixed_effects_f0_`z'_dich_size.dta", replace
export excel using "RESULTS/ins_lfabp10a/`x'/fixed_effects_f0_`z'_dich_size.xls", replace firstrow(variables)





******************************




** Merge results files with summary statistics
use "RESULTS/ins_lfabp10a/`x'/fixed_effects_f0_`y'_dich_batch.dta", clear
merge 1:1 trait model parm using "RESULTS/ins_lfabp10a/`x'/fixed_effects_f0_`y'_dich_size.dta"
drop _merge

gen order = .
replace order = 1 if strpos(parm,"1.`z'_dich")
replace order = 2 if strpos(parm,"batch")
replace order = 3 if strpos(parm,"tank")
replace order = 4 if strpos(parm,"tod")
replace order = 5 if strpos(parm,"length")
replace order = 6 if strpos(parm,"dorsal_area")
replace order = 7 if strpos(parm,"_cons")

replace trait = "intCellVol_Tot" if trait == "intVolume"
replace trait = "intCellVol_Av" if trait == "intCellVolume_Av"
replace trait = "intCellVol_Stdev" if trait == "intCellVolume_Stdev"
replace trait = "intLiver_area" if trait == "intLiverArea"
replace trait = "intLiverLip_area" if trait == "intLiverLipArea"
replace trait = "intLiverLip_n" if trait == "intLiverLipN"
replace trait = "intpigmentation" if trait == "intpig_dorsal_tail"

gen trait_order = .

replace trait_order = 1 if strpos(trait,"intNrCells")
replace trait_order = 2 if strpos(trait,"intCellVol_Tot")

replace trait_order = 5 if strpos(trait,"intCellVol_Av")
replace trait_order = 6 if strpos(trait,"intCellVol_Stdev")
replace trait_order = 7 if strpos(trait,"intCellInt_Tot")

replace trait_order = 10 if strpos(trait,"intCellInt_Av")
replace trait_order = 11 if strpos(trait,"intCellInt_Stdev")

replace trait_order = 15 if strpos(trait,"intIV")
replace trait_order = 16 if strpos(trait,"intISA")
replace trait_order = 17 if strpos(trait,"intIED")
replace trait_order = 18 if strpos(trait,"intIPAL1")
replace trait_order = 19 if strpos(trait,"intIPAL2")
replace trait_order = 20 if strpos(trait,"intIPAL3")

replace trait_order = 25 if strpos(trait,"intLiver_area")
replace trait_order = 26 if strpos(trait,"intLiverLip_n")
replace trait_order = 27 if strpos(trait,"intLiverLip_area")
replace trait_order = 28 if strpos(trait,"intLiverLipInt_Tot")
replace trait_order = 29 if strpos(trait,"intLiverLipInt_Tot_Res")

replace trait_order = 35 if strpos(trait,"intldl")
replace trait_order = 36 if strpos(trait,"inttg")
replace trait_order = 37 if strpos(trait,"inttc")
replace trait_order = 38 if strpos(trait,"intgluc")

replace trait_order = 40 if strpos(trait,"vld_area")

replace trait_order = 45 if strpos(trait,"intlength")
replace trait_order = 46 if strpos(trait,"intdorsal")
replace trait_order = 47 if strpos(trait,"intlateral")
replace trait_order = 48 if strpos(trait,"intprotein")

replace trait_order = 50 if strpos(trait,"intangle")
replace trait_order = 51 if strpos(trait,"inteye_area")
replace trait_order = 52 if strpos(trait,"intcurvature")
replace trait_order = 53 if strpos(trait,"intmajor_axis_backline")
replace trait_order = 54 if strpos(trait,"intminor_axis_backline")
replace trait_order = 55 if strpos(trait,"intpigmentation")

replace parm = "intercept" if parm == "_cons"
replace parm = "body length (in SD)" if parm == "intlength"
replace parm = "length-normalised dorsal body surface area (in SD)" if parm == "intdorsal_area"
replace parm = "10 dpf vs. 9 dpf" if parm == "10.age"
replace parm = "11 dpf vs. 9 dpf" if parm == "11.age"
replace parm = "time of imaging (per h effect)" if parm == "tod"

* make some changes to prevent work on results table

foreach var in beta_batch SE_batch LCI_batch UCI_batch t_batch P_batch N_batch      beta_size SE_size LCI_size UCI_size t_size P_size N_size {
replace `var´' = . if parm == "intercept"
}

foreach var in trait model parm {
replace `var´' = "" if parm == "intercept"
}

gen l=""

sort trait_order model order parm
drop order trait_order

format beta_batch SE_batch LCI_batch UCI_batch t_batch      beta_size SE_size LCI_size UCI_size t_size %6.3f
format P_batch P_size %12.2e
drop if model == "random" & parm == ""
drop model

order trait parm beta_batch SE_batch LCI_batch UCI_batch t_batch P_batch N_batch l     beta_size SE_size LCI_size UCI_size t_size P_size N_size

replace parm = "`z'_dich" if parm == "1.`z'_dich"
drop if strpos(parm,".batch")
drop if strpos(parm,".tank")

replace trait = substr(trait,4,.) if trait != "vld_area"


* save sumstats for Novo and for figures
outsheet using "RESULTS/ins_lfabp10a/`x'/`z'_dich_results.txt", noquote replace

}

