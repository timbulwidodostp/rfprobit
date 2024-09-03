*! Version 1.0.0 (11/22/00), G. R. Frechette  (STB-59: sg158; STB-61: sg158.1)

program define reoprob, eclass
	version 6.0
	if replay() {
		if "`e(cmd)'" ~= "reoprob" {
			error 301
		}
		Replay `0'
	}
	else	Estimate `0'
end

program define Estimate, eclass

	syntax varlist [if] [in] [, I(varname) Quadrat(int 12) /*
		*/ Level(passthru) *]

	tempvar touse x w
	marksample touse
	markout `touse' `i'

	tokenize `varlist'
	local lhs "`1'"
	macro shift 1
	local rhs "`*'"

	/* Get points and weights for Gauss-Hermite quadrature. */
	ghquadm `quadrat' `x' `w'

	/* Set up macros for ml function. */
	global S_i      "`i'"
	global S_x      "`x'"
	global S_w      "`w'"
	global S_quad   "`quadrat'"
	global S_rhs   "`rhs'"
	global S_lhs   "`lhs'"

	/* get starting values */
	tempname b0 s0
	quietly oprobit `lhs' `rhs' if `touse'
	mat `b0' = e(b)
	mat `b0' = [`b0', 0.5]
	quietly oprobit `lhs' if `touse'
	mat `s0' = e(b)
	mat `s0' = [`s0', 0.5]

	/* number of categories */
	quietly tab1 $S_lhs if `touse'
	global S_n = _result(2)

	/* create our version of `lhs' that runs from 0, ..., n-1 where
	n is the number of categories */
	tempvar dv
	rename $S_lhs `dv'
	quietly egen $S_lhs = group( `dv' ) if `touse'
	quietly replace $_lhs = $_lhs-1

	/* estimation equations */
	local meqe "($S_lhs=`rhs', nocons)"
	local start "(_cut1: $S_lhs=)"

	local i = 1
	while ( `i' < $S_n )	{
		local meqe "`meqe' /_cut`i'"
		local i = `i' + 1
	}
	local i = 2
	while ( `i' < $S_n )	{
		local start "`start' /_cut`i'"
		local i = `i' + 1
	}
	local meqe "`meqe' /rho"
	local start "`start' /rho"

	/* Sort data. */
	sort $S_i

	/* optimization */
	di in green _n "Fitting constant-only model:"
	ml model d1 reopc_ll `start' if `touse', /*
		*/ init(`s0', copy) maximize /*
		*/ search(off) /*
		*/ `options'
	di in green _n "Fitting full model:"
	ml model d1 reop_ll `meqe' if `touse', /*
		*/ continue init(`b0', copy) maximize /*
		*/ search(off) /*
		*/ `options' /*
		*/ title("Random Effects Ordered Probit")
	estimate local cmd "reoprob"
	Replay, `level'
end

program define Replay
	syntax [, Level(int $S_level)]
	ml display, level(`level')
end
