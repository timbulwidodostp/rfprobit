*! version 1.0.0  7 June 1995           sg41: STB-26
program define rfpr_ll1 /* b log_likelihood grad */ 
	version 4.0
	local b    "`1'"
	local f    "`2'"
	local g    "`3'"
	macro shift 3
	local options "FIRSTIT LASTIT FAST(string)"
	parse "`*'"

	local y    = trim("$S_mldepn")
	local doit "$S_sample"
	local i    "$S_ivar"

	tempname rho s2su beta u
	tempvar xb F p

	local rhocol = colnumb(`b',"rho:_cons")
	scalar `rho' = abs(`b'[1,`rhocol'])

	if `rho' >= 1 {
		scalar `rho' = 0.99
		di "rho >= 1, set to rho = 0.99"
	}

	matrix `b'[1,`rhocol'] = `rho'
	scalar `s2su' = sqrt(2*`rho'/(1 - `rho'))
	matrix `beta' = `b'[1,"`y':"]

	quietly {

		matrix score double `xb' = `beta' if `doit'

		gen double `F' = . in 1
		by `doit' `i': gen double `p' = cond(_n==_N,0,.) if `doit'

/* Do computation this way if only log likelihood required. */

		if "`fast'" == "0" {
			local m 1
			while `m' <= $S_quad {
				scalar `u' = `s2su'*$S_x[`m']
	
				#delimit ;
	
				by `doit' `i': replace `F' =
			    	    cond(_n==1,
					cond(`y', normprob(`xb' + `u'),
		    		      	      1 - normprob(`xb' + `u')),
					cond(`y', normprob(`xb' + `u'),
		    		      	      1 - normprob(`xb' + `u'))
					      *`F'[_n-1])
			    	    if `doit' ;

				#delimit cr
	
				replace `p' = `p' + $S_w[`m']*`F' if `doit'
	
				local m = `m' + 1
			}

			replace `F' = sum(log(`p'/sqrt(_pi))) if `doit'
			scalar  `f' = `F'[_N]

			exit
		}

/* Do computation this way if first derivatives required. */

		tempname gr
		tempvar db dr
		gen double `db' = 0
		gen double `dr' = 0

		local m 1
		while `m' <= $S_quad {
			scalar `u' = `s2su'*$S_x[`m']

			#delimit ;

			by `doit' `i': replace `F' =
			    cond(_n==1,
				cond(`y', normprob(`xb' + `u'),
		    		      1 - normprob(`xb' + `u')),
				cond(`y', normprob(`xb' + `u'),
		    		      1 - normprob(`xb' + `u'))*`F'[_n-1])
			    if `doit' ;

			replace `p' = `p' + $S_w[`m']*`F' if `doit' ;

			by `doit' `i': replace `F' = 
			    cond(`y',1,-1)*exp(-0.5*(`xb'+`u')^2)*`F'[_N]
			    /cond(`y',normprob(`xb'+`u'),1-normprob(`xb'+`u'))
			    if `doit' ;

			#delimit cr

			replace `db' = `db' + $S_w[`m']*`F'     if `doit'
			replace `dr' = `dr' + $S_w[`m']*`u'*`F' if `doit'

			local m = `m' + 1
		}

/* Compute log likelihood. */

		replace `F' = sum(log(`p'/sqrt(_pi))) if `doit'
		scalar  `f' = `F'[_N]

/* Compute first derivatives. */

		by `doit' `i': replace `p' = `p'[_N] if `doit'

		replace `db' = `db'/(sqrt(2*_pi)*`p') if `doit'

		local nvar = colsof(`beta') - 1
		if `nvar' > 0 {
			matrix `beta' = `beta'[1, 1..`nvar']
			local vars : colnames(`beta')
		}

		matrix vecaccum `g' = `db' `vars' if `doit'

		replace `dr' = sum(`dr'/`p') if `doit'

		matrix `gr' = (0)
		matrix `gr'[1,1] = `dr'[_N]/(2*sqrt(2*_pi)*`rho'*(1-`rho'))
		matrix `g' = `g' , `gr'	
	}
end
