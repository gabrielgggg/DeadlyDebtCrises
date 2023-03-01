*Main.do
*(i)    Use data from IHME to construct weekly weighted average series of 9 countries 
*(ii)   Use sovereign spread data to construct weekly weighted average spread
*(iii)  Moments using weekly data
*(iv)   Use national accounts and government debt data to generate corresponding moments

//
// change to local directory
//
cd ~/Dropbox/COVID_EM/Replication_Package/data

clear
set more off

local country "Argentina Brazil Chile Colombia Ecuador Mexico Paraguay Peru Uruguay"

****************************
**# (i) IHME data: weekly
****************************

import excel Raw_data/IHME.xls, sheet("raw_data") firstrow case(lower) clear
sort location_name date

rename location_name country


***Merge population data for these 9 countries***
merge m:1 country using "Raw_data/pop18.dta"
drop _merge


*Daily cumulative death
bys date: egen sumseir_cumulative_mean=sum(seir_cumulative_mean)
bys date: egen sumpop=sum(population)

gen avgseir_cumulative_mean=sumseir_cumulative_mean/sumpop


***Date: week of the year***
gen year=substr(date,1,4)
gen month=substr(date,6,2)
gen day = substr(date,9,2)
destring year month day, replace

gen dow = dow(mdy(month, day, year))   // generate day of week (Tuesday - February 4, 2020)
gen order=week(mdy(month, day, year))  // generate week of the year; Stata - week 1 always begins on 1 January and week 52 is always 8 or 9 days long.


***generate week order for 2020-2021 (drop 2020-02-04 - sample starts from 20200205)***
drop if date=="2020-02-04"

tostring year order, replace
forvalues i=1(1)9{ 
replace order="0`i'" if order=="`i'"
}

gen weekraw=year+order
egen weekid=group(weekraw)
sort country date
destring year order, replace


***Average mobility, mask, infection fatality***
bys country weekid: egen weeklyavg_mobility=mean(mobility_mean)
bys country weekid: egen weeklyavg_mask=mean(mask_use_mean)
bys country weekid: egen weeklyavg_infecfatality=mean(infection_fatality)


***Weighted average across countries***
bys date: asgen wavg_mobility=weeklyavg_mobility, w(population)
bys date: asgen wavg_mask=weeklyavg_mask, w(population)
bys date: asgen wavg_infecfatality=weeklyavg_infecfatality, w(population)


***Change to weely frequency***
keep country date year month day weekid wavg_mobility wavg_mask wavg_infecfatality avgseir_cumulative_mean
keep if country=="Mexico"  // Mexico data starts from 02-05-2020

bys country weekid: gen a=_n  // keep first day of each week
keep if a==1
sort date

***Generate excess death
tsset weekid
gen wavg_excessdeath=avgseir_cumulative_mean-l.avgseir_cumulative_mean
replace wavg_excessdeath=0 if wavg_excessdeath==.


preserve
keep year month day weekid
save "date_used.dta", replace   // Used for aggregate daily EMBI data to weekly frequency
restore




***************************************************
**# (i) IHME data: smooth mobility + excess death
***************************************************

*Archive: Nov 19, 2021
drop if weekid>=98


***Smooth mobility***
gen locks=max(-wavg_mobility/100,0)
replace locks=0.001 if locks==0

gen locks_orig=max(-wavg_mobility/100,0)
replace locks_orig=0.001 if locks_orig==0

*Zone 1 smooth
gen lockszone1=locks
replace lockszone1=. if weekid<9 | weekid>46

gen loglockszone1=log(lockszone1)

reg loglockszone1 c.weekid##c.weekid, noconstant
predict hatloglockszone1, xb
replace hatloglockszone1=exp(hatloglockszone1)


*Zone 2 smooth
gen lockszone2=locks
replace lockszone2=. if weekid<47 | weekid>82

gen loglockszone2=log(lockszone2)

reg loglockszone2 c.weekid##c.weekid##c.weekid##c.weekid##c.weekid, noconstant
predict hatloglockszone2, xb
replace hatloglockszone2=exp(hatloglockszone2)


*Update locks
replace locks=hatloglockszone1 if weekid>=9 & weekid<=46
replace locks=hatloglockszone2 if weekid>=47 & weekid<=82




***Smooth excess death***
gen excess=wavg_excessdeath

egen count_week=count(weekid)
gen reweight_t=weekid/count_week


forvalues i=2(1)9{ 
	gen reweight_t`i'= reweight_t^`i'
}

*OLS
mkmat reweight_t reweight_t2 reweight_t3 reweight_t4 reweight_t5 reweight_t6 reweight_t7 reweight_t8 reweight_t9, matrix(x)

mkmat excess, matrix(y)

matrix coef=inv(x'*x)*(x'*y)
matrix list coef

matrix yhat=x*coef

svmat double yhat, name(hatexcess)
replace hatexcess1=0.0 if hatexcess1<0



*Output
rename wavg_mask mask
rename wavg_infecfatality infection_fatality

keep date year month day weekid locks hatexcess1 mask infection_fatality
order date year month day weekid

save "IHME_sample_moments.dta", replace



*Export mask & smoothed series: locks + excess death
preserve
keep locks 
format locks %18.9e
export delimited data_L.txt, novarnames datafmt replace
restore

preserve
keep mask
format mask %18.9e
export delimited data_mask.txt, novarnames datafmt replace
restore

preserve
keep hatexcess1 
format hatexcess1 %18.9e
export delimited data_D.txt, novarnames datafmt replace
restore






**********************
**# (ii) Spread data
**********************
*EMBI, unit %
use "Raw_data/spread_EMBI.dta", clear


***Aggregate daily EMBI to weekly (date used in IHME - fix with nearest value if no observation on one certain day)***
merge 1:1 year month day using "date_used.dta"  // merge==2: 11/11/20; 1/1/21; 4/2/21; 12/24/21

sort year month day
foreach i of local country{
	replace `i'=`i'[_n-1] if `i'==.   // fix with nearest value if no observation on one certain day
	rename `i' embi`i'
}

drop if _merge==1  // aggregate to weekly frequency
drop Date _merge

gen date_new= mdy(month, day, year)  // date format
format date_new %d
sort date_new

order year month day weekid date_new 


***Reshape to panel***
reshape long embi, i(date_new year month day weekid) j(country) string
sort country date_new


***merge population data for these 9 countries***
merge m:1 country using "Raw_data/pop18.dta"
drop _merge


***Weighted average EMBI***
bys date_new: asgen wavg_embi=embi, w(population)
duplicates drop date_new, force


drop if weekid==98 | weekid==99

rename wavg_embi spread
save "spread_sample_moments.dta", replace








**************************
**# (iii) Moments: weekly
**************************

***Merge weekly sample***
use "IHME_sample_moments.dta", clear

merge 1:1 weekid using "spread_sample_moments.dta"
drop _merge


*Base: spread
quietly su spread if weekid==1, detail
local sp_1= r(mean)


*Sample: since April 1, 2020
keep if weekid>=9


*Average case fatality rate
egen mean_case_fataility = mean(infection_fatality)
replace mean_case_fataility=round(100*mean_case_fataility, 0.01)


*Moments: cumulative fatalities
sort weekid
gen cum_excess=hatexcess1[1]
replace cum_excess=hatexcess1[_n]+ cum_excess[_n-1] if _n>1


*Moments: cumulative fatalities 2020/2021 (per 100)
quietly su cum_excess if year == 2020, detail
local mean_cumexcess2020 = r(mean)

quietly su cum_excess if year == 2021, detail
local mean_cumexcess2021 = r(mean)

gen mean_cumexcess_2020 = round(100*`mean_cumexcess2020', 0.01)  // mean 2020
gen mean_cumexcess__2021 = round(100*`mean_cumexcess2021', 0.01)  // mean 2021


*Moments: social distancing
quietly su locks if year == 2020, detail
local mean_lock_2020 = r(mean)
local max_lock_2020 = r(max)

quietly su locks if year == 2021, detail
local mean_lock_2021 = r(mean)
local max_lock_2021 = r(max)

gen max_lock_2020  = round(`max_lock_2020', 0.01)       // max 2020
gen mean_lock_2020 = round(`mean_lock_2020', 0.01)     // mean 2020
gen max_lock_2021  = round(`max_lock_2021', 0.01)     // max 2021
gen mean_lock_2021 = round(`mean_lock_2021', 0.01)     // mean 2021


*Moments: Spread 
quietly su spread if year == 2020, detail
local mean_sp_2020 = r(mean)
local max_sp_2020 = r(max)

gen sp_peak_rise_2020 = round(`max_sp_2020'-`sp_1', 0.1)
gen sp_mean_rise_2020 = round(`mean_sp_2020'-`sp_1', 0.1)



***Output***
collapse (mean) mean_case_fataility mean_cumexcess_2020 mean_cumexcess__2021 max_lock_2020 mean_lock_2020 max_lock_2021 mean_lock_2021 sp_peak_rise_2020 sp_mean_rise_2020
order mean_case_fataility mean_cumexcess_2020 mean_cumexcess__2021 max_lock_2020 mean_lock_2020 max_lock_2021 mean_lock_2021 sp_peak_rise_2020 sp_mean_rise_2020

xpose, clear varname format
rename _varname economic_moments
rename v1 moments
order economic_moments moments


*Format output
gen order=_n
replace order=order+3 if economic_moments=="sp_peak_rise_2020" | economic_moments=="sp_mean_rise_2020"

count
set obs `=`r(N)'+3'
replace order=8 if _n==10
replace order=9 if _n==11
replace order=10 if _n==12

sort order
order order


export excel using  "Output.xls",  sheet("Table 1", modify) firstrow(variables) 






**********************************************
**# (iv) Government debt/GDP + consumption
**********************************************

use "Raw_data/NA_govtdebt.dta", clear


***Merge population data for these 9 countries***
merge m:1 country using "Raw_data/pop18.dta"
drop _merge


***Convert debt (now USD mn) into local currency: LC mn***
*Exchange_rate "Official Rate: End of Period: National Currency per USD"
gen govtdebt_LC=govtdebt*exchange_rate

*Real debt = nominal debt/ (nominal GDP sa/real GDP sa)
gen realgovtdebt_LC=100*govtdebt_LC/GDPdeflator_sa2019


***Debt-to-GDP(GDP2019) ratio (LC mn/LC mn)***
gen govtdebt_GDP2019 = 100*realgovtdebt_LC/annualGDP_nom2019


***Consumption decline 2020***
*Consumption: private + public
gen consumption=C+G

*Annual consumption in 2019 and 2020
gen annual_consumption2019=.
gen annual_consumption2020=.

foreach i of local country{
	quietly su consumption if year == 2019 & country=="`i'", detail
	local meanCG2019 = r(mean)
	replace annual_consumption2019=`meanCG2019' if country=="`i'"
	
	quietly su consumption if year == 2020 & country=="`i'", detail
	local meanCG2020 = r(mean)
	replace annual_consumption2020=`meanCG2020' if country=="`i'"
}


***Weighted average consumption***
preserve

sort country quarterdate
gen growthCG2020=100*(ln(annual_consumption2020)-ln(annual_consumption2019))

keep country growthCG2020 population
duplicates drop country, force
asgen Consumption_decline_2020=growthCG2020, w(population)

local Consumption_decline_2020 = Consumption_decline_2020[1]
restore

gen Consumption_decline_2020=round(`Consumption_decline_2020', 0.1)


***Weighted average debt-to-GDP ratio***
sort country quarterdate
bys quarterdate: asgen wavggovtdebt_GDP2019=govtdebt_GDP2019, w(population)

duplicates drop quarterdate, force
drop if year==.
keep year quarterdate Consumption_decline_2020 wavggovtdebt_GDP2019 


***Government debt to output 2019***
quietly su wavggovtdebt_GDP2019 if year == 2019, detail
local mean2019 = r(mean)
gen Govtdebt_output2019=round(`mean2019',1)


***Government debt increase 2020***
quietly su wavggovtdebt_GDP2019 if quarterdate==tq(2019,4), detail
local end2019 = r(mean)
gen end2019=`end2019'

quietly su wavggovtdebt_GDP2019 if year==2020, detail
local max2020 = r(max)
gen max2020=`max2020'

gen Govtdebt_increase_2020=round(max2020-end2019, 1)


***Output***
collapse (mean) Govtdebt_output2019 Govtdebt_increase_2020 Consumption_decline_2020
order Govtdebt_output2019 Govtdebt_increase_2020 Consumption_decline_2020

xpose, clear varname format
rename _varname economic_moments
rename v1 moments
order economic_moments moments

export excel using Output.xls, cell(B9) sheet("Table 1", modify)



**********************
**# House keeping
**********************
capture erase "IHME_sample_moments.dta"
capture erase "spread_sample_moments.dta"
capture erase "date_used.dta"
