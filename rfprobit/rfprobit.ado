*! version 1.1.0  20jun1995     sg41: STB-26
program define rfprobit
	version 4.0
	local options "Level(integer $S_level)"
	if substr("`1'",1,1)=="," | "`*'"=="" { 
		if "$S_E_cmd"~="rfprobit" { 
			error 301
		}
		parse "`*'"
	}
	else { 
		local varlist "req ex"
		local if "opt"
		local in "opt"
		local options /*
		*/ "`options' I(string) Quadrat(integer 6) noCHIsq noLOg *"
		parse "`*'"
		parse "`varlist'", parse(" ")
		local y "`1'"
		macro shift
		xt_iis `i'
		local i "$S_1"
		tempvar doit x w
		mark `doit' `if' `in'
		markout `doit' `y' `*' `i'

	/* Check to see if outcome varies. */

		quietly count if `doit'
		local n = _result(1)
		quietly count if `y'==0 & `doit'
		local n0 = _result(1)
		if `n0'==0 | `n0'==`n' {
			di _n in blu "outcome does not vary"
			exit
		}

	/* Sort data. */

		sort `doit' `i'

	/* Get points and weights for Gaussian-Hermite quadrature. */
		
		ghquad double(`x' `w'), n(`quadrat')

	/* Set up macros for ml function. */

		global S_sample "`doit'"
		global S_ivar   "`i'"
		global S_x      "`x'"
		global S_w      "`w'"
		global S_quad   "`quadrat'"

	/* Fit constant-only model. */

		if "`chisq'"=="" & "`*'"~="" {
			rfpr_ml `y', title("Constant-only model") /*
			*/	     `log' `options'
			local lf0 "lf0($S_1)"
			local rho "rho($S_2)"
		}

	/* Fit full model. */

		rfpr_ml `y' `*', title("Full model") post `lf0' `rho' /*
		*/	         `log' `options'
	}

/* Display results. */

	ml mlout rfprobit, level(`level')

/* Compute LR test for rho = 0. */

	local chi = 2*($S_E_ll - $S_E_rho0)

	#delimit ;
	di in gr "LR test of rho = 0:   chi2(" in ye "1" in gr ")     = "
	   in ye %7.2f `chi' _n
	   in gr "                      Prob > chi2 = " in ye %7.4f
	   chiprob(1,`chi') ;
	#delimit cr
end


program define rfpr_ml /* y x, title(string) POST noLOg ml_options */ 
	version 4.0
	local varlist "req ex"
	local options "TITLE(string) POST LF0(string) RHO(string) noLOg *"
	parse "`*'"
	parse "`varlist'", parse(" ")
	local y "`1'"
	macro shift
	local doit "$S_sample"
	if "`lf0'"~="" { local lf0 "lf0(`lf0')" }

	tempname llrho0 b0 b1 lllast b ll V
	tempvar  mldoit

/* Get initial values. */
		
	quietly probit `y' `*' if `doit'
	scalar `llrho0' = _result(2)
	if "`log'"=="" {
		di _n in gr "`title'"
		di in gr "rho =" in ye %4.1f 0 /*
		*/ in gr "     Log Likelihood = " in ye `llrho0'
	}
	matrix `b0' = get(_b)
	matrix coleq `b0' = `y'
	matrix `b1' = (0)
	matrix colnames `b1' = rho:_cons
	matrix `b0' = `b0' , `b1'
	local rcol = colnumb(`b0',"rho:_cons")

/* Search for good starting value for rho if not supplied. */

	if "`rho'"=="" {
		global S_mldepn "`y'"
		scalar `lllast' = `llrho0'
		local rho    0.05
		local rhotry 0.1
		while `rhotry' < 0.91 {
			matrix `b0'[1,`rcol'] = `rhotry'
			rfpr_ll1 `b0' `ll' `b' , fast(0)
			if "`log'"=="" {
				di in gr "rho =" in ye %4.1f `rhotry' /*
				*/ in gr "     Log Likelihood ~ " in ye `ll'
			}
			if `ll' < `lllast' { /* exit loop */
				local rhotry 1
			}
			else {
				scalar `lllast' = `ll'
				local rho `rhotry'
				local rhotry = `rhotry' + 0.1
			}
		}
	}

	matrix `b0'[1,`rcol'] = `rho'

/* Set up ml commands. */

	ml begin
	ml function rfpr_ll1
	ml method deriv1
	eq `y': `y' `*'
	eq rho:
	ml model `b' = `y' rho, depv(10) from(`b0')
	ml sample `mldoit' if `doit', noauto

	if "`log'"=="" { noisily ml max `ll' `V', `options' }
	else		 quietly ml max `ll' `V', `options'

	global S_1 = `ll'
	global S_2 = `b'[1,`rcol']

	if "`post'"~="" {
		ml post rfprobit, `lf0' pr2 title("Random-Effects Probit")
	}

	global S_E_rho0 = `llrho0'
end

/*
	Routines that compute weights and points for Gaussian-Hermite
	quadrature follow:
*/

* version 1.0.1  29jun1995
program define ghquad 
	version 4.0
	local varlist "req new min(2) max(2)"
	local options "N(integer 10)"
	parse "`*'"
	parse "`varlist'", parse(" ")
	local x "`1'"
	local w "`2'"
	if `n' + 2 > _N  {
		di in red  /*
		*/ "`n' + 2 observations needed to compute quadrature points"
		exit 2001
	}
	tempname xx ww
	local i 1
	local m = int((`n' + 1)/2)
	while `i' <= `m' {
		if `i' == 1 {
			scalar `xx' = sqrt(2*`n'+1)-1.85575*(2*`n'+1)^(-1/6)
		}
		else if `i' == 2 { scalar `xx' = `xx'-1.14*`n'^0.426/`xx' }
		else if `i' == 3 { scalar `xx' = 1.86*`xx'-0.86*`x'[1] }
		else if `i' == 4 { scalar `xx' = 1.91*`xx'-0.91*`x'[2] }
		else { scalar `xx' = 2*`xx'-`x'[`i'-2] }
		hermite `n' `xx' `ww'
		qui replace `x' = `xx' in `i'
		qui replace `w' = `ww' in `i'
		local i = `i' + 1
	}
	if mod(`n', 2) == 1 { qui replace `x' = 0 in `m' }
	qui replace `x' = -`x'[`n'+1-_n] in `i'/`n'
	qui replace `w' =  `w'[`n'+1-_n] in `i'/`n'
end


program define hermite  /* integer n, scalar x, scalar w */
	version 4.0
	local n "`1'"
	local x "`2'"
	local w "`3'"
	local last = `n' + 2
	tempvar p
	tempname i
	qui gen double `p' = . 
	scalar `i' = 1
	while `i' <= 10 {
		qui replace `p' = 0 in 1
		qui replace `p' = _pi^(-0.25) in 2
		qui replace `p' = `x'*sqrt(2/(_n-2))*`p'[_n-1] /*
		*/	- sqrt((_n-3)/(_n-2))*`p'[_n-2] in 3/`last'
		scalar `w' = sqrt(2*`n')*`p'[`last'-1]
		scalar `x' = `x' - `p'[`last']/`w'
		if abs(`p'[`last']/`w') < 3e-14 {
			scalar `w' = 2/(`w'*`w')
			exit
		}
		scalar `i' = `i' + 1
	}
	di in red "hermite did not converge"
	exit 499
end
	
