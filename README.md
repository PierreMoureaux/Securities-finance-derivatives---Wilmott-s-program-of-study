Paul Wilmott’s program of study

I – Introduction
In his famous book http://iman.sn/bibliotek/livres/filieres/banque-finance-assurance/pdfs/paul-wilmott-on-quantitative-finance.pdf, Paul Wilmott highlighted several program of studies about numerical methods (especially finite difference and monte-carlo simulation). The purpose of this git-hub stream will be to follow this program and show numerical results

II – Programs of study
Programs of study are following ones:
-	1/Finite difference
   - Explicit method/European calls, puts and binaries: To get started you should learn the explicit method as applied to the Black–Scholes equation for a European option.  
   - Explicit method/American calls, puts and binaries  
   - Crank–Nicolson/European calls, puts and binaries
   - Crank–Nicolson/American calls, puts and binaries
   - Explicit method/path-dependent options 
   - Interest rate products
   - Two-factor explicit
   - Two-factor implicit
-	2/Monte Carlo
- European calls, puts and binaries on a single equity
- Path-dependent option on a single equity: Price a barrier, Asian, lookback etc.
- Options on many stocks: Price a multi-asset contract by simulating correlated random walks
- Interest rate derivatives, spot rate model
- HJM model
- BGM model: A discrete version of HJM.
- Bonus : Generalized RFR BGM model

III – Disclaimers
-	For each stream, a summary documentation will be attached (merely mathematical oriented)
-	Targets are Python and C++ codes, and especially numerical results
-	Sources will be documented (of course Paul Wilmott’s book but also all other sources)
-	Python and C++ codes will be assumed to be improved, but I’ll try to do my best in order to be in line with modern python and C++ concepts, without “old way” programming
