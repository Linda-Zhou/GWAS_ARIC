insheet using pc.txt, clear
rename sampleid iid
save pc.dta, replace

insheet using logNBL1_SeqId_2944_66.visit3.W.FAST.runfile.iid.txt, clear

drop if log=="NA"
drop if pc1=="NA"
drop pc*

merge 1:1 iid using pc.dta, keep(match)
foreach var of varlist * {
  destring `var', replace
 }

regress lognbl1_seqid_2944_66  age gender pc* forsyth minneapolis washington
predict res_f, residuals

egen rank=rank(res_f)
su res_f, meanonly
gen invres_f=invnormal((rank-0.5)/r(N))


regress lognbl1_seqid_2944_66  age  pc* forsyth minneapolis washington
predict res, residuals

egen rank1=rank(res)
su res, meanonly
gen invres=invnormal((rank1-0.5)/r(N))


outsheet using residuals.txt, replace
