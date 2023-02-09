Paul Wilmott’s program of study

I – Introduction
In his famous book http://iman.sn/bibliotek/livres/filieres/banque-finance-assurance/pdfs/paul-wilmott-on-quantitative-finance.pdf, Paul Wilmott highlighted several program of studies about numerical methods (especially finite difference and monte-carlo simulation). The purpose of this git-hub stream will be to follow this program and show numerical results

II – Programs of study
Programs of study are following ones:
-	1/Finite difference
   - Explicit method/European calls, puts and binaries: To get started you should learn the explicit method as applied to the Black–Scholes equation for a European option. This is very easy to program and you won’t make many mistakes. 
   - Explicit method/American calls, puts and binaries: Not much harder is the application of the explicit method to American options.  
   - Crank–Nicolson/European calls, puts and binaries: Once you’ve got the explicit method under your belt you should learn the Crank–Nicolson implicit method. This is harder to program, but you will get a better accuracy.
   - Crank–Nicolson/American calls, puts and binaries: There’s not much more effort involved in pricing American-style options than in the pricing of European-style options.  
   - Explicit method/path-dependent options: By now you’ll be quite sophisticated and it’s time to price a path-dependent contract. Start with an Asian option with discrete sampling, and then try a continuously sampled Asian. Finally, try your hand at lookbacks.  
   - Interest rate products: Repeat the above program for non-path-dependent and then path dependent interest rate products. First price caps and floors and then go on to the index amortizing rate swap. 
   - Two-factor explicit: To get started on two-factor problems price a convertible bond using an explicit method, with both the stock and the spot interest rate being stochastic.  
   - Two-factor implicit: The final stage, using methods described in this book, is to implement the implicit two-factor method as applied to the convertible bond.
-	2/Monte Carlo
- European calls, puts and binaries on a single equity: Simulate a single stock path, the payoff for an option, or even a portfolio of options, calculate the expected payoff and present value to price the contract.
- Path-dependent option on a single equity: Price a barrier, Asian, lookback etc.
- Options on many stocks: Price a multi-asset contract by simulating correlated random walks. You’ll see how time taken varies with number of dimensions.
- Interest rate derivatives, spot rate model: This is not that much harder than equities. Just remember to present value along each realized path of rates before taking the expectation across all paths.
- HJM model: Slightly more ambitious is the HJM interest rate model. Use a single factor, then two factors etc.
- BGM model: A discrete version of HJM.
- Bonus : Generalized RFR BGM model

III – Disclaimers
-	For each stream, a summary documentation will be attached (merely mathematical oriented)
-	Targets are Python and C++ codes, and especially numerical results
-	Sources will be documented (of course Paul Wilmott’s book but also all other sources)
-	Python and C++ codes will be assumed to be improved, but I’ll try to do my best in order to be in line with modern python and C++ concepts, without “old way” programming
