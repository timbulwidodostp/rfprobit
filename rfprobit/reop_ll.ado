*! Version 1.0.1 November 22 2000, by Guillaume R. Frechette  (STB-59: sg158)

program define reop_ll
	version 6.0
	args todo b lnf g

	tempvar theta1 F F1 F2 p db db1 db2 dr C
	tempname rho s2su u x w gr
	local nm1 = $S_n-1
	local np1 = $S_n+1
	local i = 1
	while `i' < $S_n	{
		tempname _cut`i'
		tempname g_cut`i'
		tempname g_cut1`i'
		tempname g_cut2`i'
		local i = `i'+1
	}
	mleval `theta1' = `b', eq(1)
	local h 0
	local i 1
	local j 2
	while `i' < $S_n	{
		mleval `_cut`i'' = `b', eq(`j') scalar
		if `h'>0 {
			if `_cut`i''<=`_cut`h'' {
				scalar `_cut`h''=`_cut`i''-0.01
				/* the above correction should not be needed
				in general, but it might preclude `bad'
				estimates if, for instance, one performs a
				random search */
			}
		} 
		local h = `h'+1
		local i = `i'+1
		local j = `j'+1
	}
	mleval `rho' = `b', eq(`np1') scalar

	if `rho' >= 1 {
		scalar `rho' = 0.99
		di "rho >= 1, set to rho = 0.99"
	}

	scalar `s2su' = sqrt(2*`rho'/(1 - `rho'))

	quietly {

		gen double `F' = . in 1
		gen double `F1' = . in 1
		gen double `F2' = . in 1
		by $S_i: gen double `p' = cond(_n==_N,0,.)
		gen double `db' = 0
		gen double `db1' = 0
		gen double `db2' = 0
		gen double `dr' = 0
		gen double `C' = .

		local m 1
		while `m' <= $S_quad {
			scalar `x' = $S_x[1,`m']
			scalar `w' = $S_w[1,`m']
			scalar `u' = `s2su'*`x'

			local condf "cond($S_lhs==0, normprob(`_cut1'-`theta1'-`u')"
			local i 1
			local j 2
			while `j' < $S_n {
				local condf "`condf', cond($S_lhs==`i', normprob(`_cut`j''-`theta1'-`u')-normprob(`_cut`i''-`theta1'-`u')"
				local i = `i'+1
				local j = `j'+1
			}
			local condf "`condf', 1-normprob(`_cut`nm1''-`theta1'-`u'))"
			local i 2
			while `i' < $S_n {
				local i = `i'+1
				local condf "`condf')"
			}
			replace `C' = `condf'
			replace `C' = 0.00000001 if `C' == 0

			by $S_i: replace `F' = /*
				*/ cond(_n==1,`C',`C'*`F'[_n-1])

			replace `p' = `p' + `w'*`F'

			local condg1 "cond($S_lhs==0, 0"
			local condg2 "cond($S_lhs==0, -exp(-0.5*(`_cut1'-`theta1'-`u')^2)"
			local h 1
			local i 2
			while `i' < $S_n {
				local condg1 "`condg1', cond($S_lhs==`h', exp(-0.5*(`_cut`h''-`theta1'-`u')^2)"
				local condg2 "`condg2', cond($S_lhs==`h', -exp(-0.5*(`_cut`i''-`theta1'-`u')^2)"
				local h = `h'+1
				local i = `i'+1
			}
			local condg1 "`condg1', exp(-0.5*(`_cut`nm1''-`theta1'-`u')^2))"
			local condg2 "`condg2', 0)"
			local i 2
			while `i' < $S_n {
				local i = `i'+1
				local condg1 "`condg1')"
				local condg2 "`condg2')"
			}

			replace `F1' = `F'
			replace `F2' = `F'
			by $S_i: replace `F1' = `condg1'*`F1'[_N]/`C'
			by $S_i: replace `F2' = `condg2'*`F2'[_N]/`C'
			by $S_i: replace `F' = (`condg1'+`condg2')*`F'[_N]/`C'

			replace `db' = `db' + `w'*`F'
			replace `db1' = `db1' + `w'*`F1'
			replace `db2' = `db2' + `w'*`F2'
			replace `dr' = `dr' + `w'*`u'*`F'

			local m = `m' + 1
		}

		tempname lp grps
		gen double `lp' = ln(`p'/sqrt(_pi))
		by $S_i: gen byte `grps'=_n==_N
		sum `grps' if $ML_samp
		local N = r(sum)
		sum `lp' if $ML_samp, meanonly
		if r(N) !=`N' {
			scalar `lnf' = .
			exit
		}
		scalar `lnf' = r(sum)

		if `todo'==0|`lnf'==. { exit }

		by $S_i: replace `p' = `p'[_N]

		replace `db' = `db'/(sqrt(2*_pi)*`p')
		replace `db1' = -`db1'/(sqrt(2*_pi)*`p')
		replace `db2' = -`db2'/(sqrt(2*_pi)*`p')

		matrix vecaccum `g' = `db' $S_rhs, nocons
		local i 0
		local j 1
		while `j' < $S_n	{
		capture {
				matrix vecaccum `g_cut1`j'' = `db1' if $S_lhs == `j'
				matrix vecaccum `g_cut2`j'' = `db2' if $S_lhs == `i'
			}
			if _rc {
				scalar `lnf' = .
				exit
			}
			matrix `g_cut`j'' = `g_cut1`j''+`g_cut2`j''
			matrix `g' = `g', `g_cut`j''
			local i = `i'+1
			local j = `j'+1
		}

		replace `dr' = sum(`dr'/`p')

		scalar `gr' = `dr'[_N]/(2*sqrt(2*_pi)*`rho'*(1-`rho'))
		matrix `g' = `g', `gr'
	}
end
