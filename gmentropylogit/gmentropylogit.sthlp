{smcl}
{* *! version 1.0.0  27nov20013}{...}
{cmd:help gmentropylogit}{right: ({browse "http://www.stata-journal.com/article.html?article=up0050":SJ16-1: st0390_1})}
{hline}

{title:Title}

{p2colset 5 23 25 2}{...}
{p2col:{hi: gmentropylogit} {hline 2}}Generalized maximum entropy discrete choice model{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 11 2}
{opt gmentropylogit:}
{depvar} [{indepvars}] {ifin}
[{cmd:,} {opt m:fx} {opt gen:erate(varname)} {opt p:riors(varlist numeric max=1)}]


{title:Description}

{pstd}
Given finite samples, {cmd:gmentropylogit} is more efficient than
its maximum entropy and maximum likelihood counterparts because it
incorporates noise terms in its results.  It also performs better than
its maximum likelihood counterparts when working with small sample
sizes.  {cmd:gmentropylogit} models the probability of a positive
outcome given a set of regressors.  {it:depvar} is equal to zero or one.
{it:depvar}={cmd:0} indicates a negative outcome, and
{it:depvar}={cmd:1} indicates a positive outcome.


{title:Options}

{phang}
{cmd:mfx} displays, instead of the coefficients, the impact of each x on the
probability of a positive outcome.  {cmd:mfx} considers dummy variables and
provides their estimates accordingly.

{phang}
{cmd:generate(}{it:varname}{cmd:)} creates a new variable with a user-selected
name that will contain the predicted probability of the fitted model.

{phang}
{cmd: priors(}{it:varlist numeric max=1}{cmd:)} Optional, specifies a variable with
priors for each observation.


{title:Example}

{phang}{cmd:. webuse lbw}{p_end}
{phang}{cmd:. gmentropylogit low age smoke ht, mfx}{p_end}


{marker Authors}{...}
{title:Authors}

{pstd}Paul Corral{p_end}
{pstd}American University{p_end}
{pstd}Washington, DC{p_end}
{pstd}paulcorral@gmail.com{p_end}

{pstd}Mungo Terbish{p_end}
{pstd}American University{p_end}
{pstd}Washington, DC{p_end}
{pstd}mungunsuvd@gmail.com{p_end}


{marker also_see}{...}
{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 16, number 1: {browse "http://www.stata-journal.com/article.html?article=up0050":st0390_1},{break}
                    {it:Stata Journal}, volume 15, number 2: {browse "http://www.stata-journal.com/article.html?article=st0390":st0390}
{p_end}
